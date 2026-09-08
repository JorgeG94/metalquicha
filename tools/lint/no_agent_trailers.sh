#!/usr/bin/env bash
# commit-msg hook: keep AI-agent session trailers out of commit messages.
#
# Agents that write commits here are told to append a link back to the
# session that produced them.  That link is useful in a chat transcript
# and worthless in `git log`: it points at a conversation nobody outside
# the session can open, it never resolves again once the session is
# archived, and it is noise in every `git log`, `git blame` and release
# note from here on.  The repository is the durable artefact; a session
# URL is not.
#
# Enforced rather than written down because the instruction to add the
# trailer lives in the agent's own configuration, which this repository
# does not control and cannot see.  A rule in AGENTS.md competes with
# that configuration and loses; a hook that rejects the commit does not.
#
# Two modes, one set of patterns:
#
#     tools/lint/no_agent_trailers.sh .git/COMMIT_EDITMSG
#     tools/lint/no_agent_trailers.sh --range origin/main..HEAD
#
# The first is the commit-msg hook, and is what pre-commit runs.  The
# second re-checks history in CI, because a commit-msg hook only fires
# for someone who has run `pre-commit install`, and `--no-verify` skips
# it outright.  Neither is a substitute for the other.
set -euo pipefail

# Anchored where it can be, so a commit that *discusses* the convention --
# this one, for instance -- is not caught by its own rule.  What is banned
# is a trailer or a bare session link, not the words.
PATTERNS=(
    '^[[:space:]]*Claude-Session[[:space:]]*:'
    '^[[:space:]]*Session-Link[[:space:]]*:'
    'https?://[a-z.]*claude\.ai/code/session'
    'Generated with \[Claude Code\]'
)

usage() {
    cat >&2 <<'EOF'
usage: no_agent_trailers.sh <commit-msg-file>
       no_agent_trailers.sh --range <git-rev-range>
EOF
    exit 2
}

# Prints the offending lines of a message, empty if it is clean.  One
# line can match several patterns -- a trailer that is also a URL matches
# two -- so the hits are uniqued on line number rather than reported once
# per pattern.
offending() {
    local message=$1 pattern hits=""
    for pattern in "${PATTERNS[@]}"; do
        hits+=$(printf '%s\n' "$message" | grep -nE "$pattern" || true)$'\n'
    done
    printf '%s' "$hits" | grep -v '^$' | sort -n -u -t: -k1,1 || true
}

explain() {
    cat >&2 <<'EOF'

ERROR: this commit message carries an AI-agent session trailer.

Session links belong in the conversation, not in git history: they are
unopenable by anyone else, they rot when the session is archived, and
they outlive their usefulness in every later `git log`.

Drop the offending line and commit again.  To fix a commit already made:

    git commit --amend            # the most recent one
    git rebase --exec 'git commit --amend --no-edit' <base>

Co-authorship trailers are fine and are not what this checks.
EOF
}

status=0

case "${1:-}" in
    "") usage ;;
    --range)
        [[ $# -eq 2 ]] || usage
        # %B is the raw body, so a trailer split across a wrapped line is
        # still seen the way git stored it.
        while read -r sha; do
            [[ -n "$sha" ]] || continue
            hits=$(offending "$(git log -1 --format=%B "$sha")")
            if [[ -n "$hits" ]]; then
                printf '%s: %s\n' "${sha:0:12}" \
                    "$(git log -1 --format=%s "$sha")" >&2
                printf '  %s\n' "$hits" >&2
                status=1
            fi
        done < <(git rev-list "$2")
        ;;
    -*) usage ;;
    *)
        [[ -f "$1" ]] || { echo "no such file: $1" >&2; exit 2; }
        # Comment lines are stripped by git before the message is stored,
        # so a `#` line mentioning the trailer must not fail the commit.
        hits=$(offending "$(grep -v '^#' "$1" || true)")
        if [[ -n "$hits" ]]; then
            printf '%s\n' "$hits" >&2
            status=1
        fi
        ;;
esac

[[ $status -eq 0 ]] || explain
exit $status
