"""Strict validation helpers for MQC input parsing."""

from __future__ import annotations
from typing import Any, Dict


def die(msg: str) -> None:
    """Raise a ValueError with the given message."""
    raise ValueError(msg)


def require_only_keys(d: Dict[str, Any], allowed: set[str], ctx: str) -> None:
    """Ensure dictionary only contains allowed keys."""
    extra = set(d.keys()) - allowed
    if extra:
        raise ValueError(f"{ctx}: unknown keys: {sorted(extra)}")


def req_type(x: Any, t, ctx: str) -> Any:
    """Require that x is of type t, raising ValueError if not."""
    if not isinstance(x, t):
        die(f"{ctx}: expected {t.__name__}, got {type(x).__name__}")
    return x


def opt_type(x: Any, t, ctx: str) -> Any:
    """Check type if x is not None."""
    if x is None:
        return None
    if not isinstance(x, t):
        die(f"{ctx}: expected {t.__name__}, got {type(x).__name__}")
    return x


def check_index_0based(i: Any, natom: int, ctx: str) -> int:
    """Validate that i is a valid 0-based atom index."""
    if not isinstance(i, int):
        die(f"{ctx}: expected int, got {type(i).__name__}")
    if i < 0 or i >= natom:
        die(f"{ctx}: index out of range {i} for natom={natom} (0-based expected)")
    return i
