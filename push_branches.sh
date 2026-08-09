#!/usr/bin/env bash
#
# Push the PR stack, in dependency order.
#
# Each branch is a contiguous prefix of feat/json-input-followup, so PR n
# targets PR n-1 rather than main. They must go up in this order: GitHub will
# not let you set a base branch that does not exist yet.
#
#   ./push_branches.sh            show what would happen, push nothing
#   ./push_branches.sh --push     actually push
#   ./push_branches.sh --push --only feat/cpu-backend    just one
#
# Nothing here force-pushes. A branch whose remote has moved is reported and
# skipped, because the only reasons for that are someone else's work or a
# rewrite, and neither should be resolved by a script.

set -uo pipefail

REPO="${REPO:-/home/jorge/dev/metalquicha}"
REMOTE="${REMOTE:-origin}"
BASE="${BASE:-main}"

# In order. The second field is the base a PR from this branch should target.
STACK=(
  "fix/gmbe-thread-pin|main|one commit, independent of the rest"
  "feat/json-input|main|breaking 0.2.0: retire the .mqc format"
  "feat/ecp-reader|feat/json-input|read ECPs from BSE files"
  "feat/python-interface|feat/ecp-reader|C API, MPI session, Python package"
  "fix/scf-convergence|feat/python-interface|SCF convergence, GMBE pin, allow_crap_scf"
  "feat/checkpoint-restart|fix/scf-convergence|fingerprint and text checkpoints"
  "feat/hdf5-checkpoints|feat/checkpoint-restart|HDF5, derivative checkpoints"
  "feat/properties|feat/hdf5-checkpoints|HOMO/LUMO, n-mer multiplicity"
  "feat/cpu-backend|feat/properties|libcint, CPU RHF, the contraction fix"
  "feat/cpu-fragmented-hf|feat/cpu-backend|fragmented HF, checked against PySCF"
  "feat/density-fitting|feat/cpu-fragmented-hf|density-fitted J and K"
  "feat/libxc|feat/density-fitting|libxc for functionals"
)

DO_PUSH=0
ONLY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --push) DO_PUSH=1 ;;
    --only) ONLY="${2:-}"; shift ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done

cd "$REPO" || exit 1

if [[ -n "$(git status --porcelain)" ]]; then
  echo "working tree is dirty -- commit or stash first" >&2
  exit 1
fi

echo "remote: $REMOTE  ($(git remote get-url "$REMOTE"))"
[[ $DO_PUSH -eq 1 ]] || echo "DRY RUN -- nothing will be pushed. Add --push when the list looks right."
echo

pushed=(); skipped=(); blocked=()

for entry in "${STACK[@]}"; do
  IFS='|' read -r branch base note <<<"$entry"
  [[ -n "$ONLY" && "$ONLY" != "$branch" ]] && continue

  if ! git show-ref --verify --quiet "refs/heads/$branch"; then
    echo "  MISSING   $branch  (no such local branch)"
    blocked+=("$branch")
    continue
  fi

  local_sha=$(git rev-parse "$branch")
  remote_sha=$(git ls-remote --heads "$REMOTE" "$branch" 2>/dev/null | awk '{print $1}')
  n=$(git rev-list --count "$BASE..$branch")

  if [[ -z "$remote_sha" ]]; then
    state="new"
  elif [[ "$remote_sha" == "$local_sha" ]]; then
    echo "  UP TO DATE $branch  (${n} commits)  -> PR base: $base"
    skipped+=("$branch")
    continue
  elif git merge-base --is-ancestor "$remote_sha" "$local_sha"; then
    state="fast-forward"
  else
    # The remote has commits this branch does not. Pushing would need --force,
    # which would discard them.
    echo "  DIVERGED  $branch  -- remote has work not in the local branch; skipping"
    blocked+=("$branch")
    continue
  fi

  if [[ $DO_PUSH -eq 1 ]]; then
    if git push "$REMOTE" "$branch"; then
      echo "  PUSHED    $branch  ($state, ${n} commits)  -> PR base: $base"
      pushed+=("$branch")
    else
      echo "  FAILED    $branch"
      blocked+=("$branch")
    fi
  else
    printf "  would push %-24s %-13s %3s commits  -> PR base: %s\n" \
      "$branch" "($state)" "$n" "$base"
    printf "             %s\n" "$note"
  fi
done

echo
echo "pushed ${#pushed[@]}, already current ${#skipped[@]}, skipped ${#blocked[@]}"

if [[ $DO_PUSH -eq 1 && ${#pushed[@]} -gt 0 ]]; then
  cat <<'NOTE'

Open the PRs in the same order -- each targets the one before it:

NOTE
  for entry in "${STACK[@]}"; do
    IFS='|' read -r branch base note <<<"$entry"
    [[ -n "$ONLY" && "$ONLY" != "$branch" ]] && continue
    printf '  gh pr create --base %-24s --head %-24s --fill\n' "$base" "$branch"
  done
fi
