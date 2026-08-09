"""Loading libmqc and declaring its C entry points.

Kept apart from the classes so that the ctypes -- which is the part that is
tedious and the part that segfaults when it is wrong -- is in one place and
readable as a list.

Note the library is loaded with RTLD_GLOBAL. It brings its own copy of the MPI
runtime, and a second one loaded privately would have the ranks talking past
each other.
"""

import ctypes
import os
import sys

_ENV = "MQC_LIBRARY"

_c_int = ctypes.c_int
_c_long = ctypes.c_int64
_c_double = ctypes.c_double
_c_ptr = ctypes.c_void_p
_c_str = ctypes.c_char_p


def _candidates():
    """Where libmqc might be, most explicit first."""
    if os.environ.get(_ENV):
        yield os.environ[_ENV]
        return
    name = "libmqc.dylib" if sys.platform == "darwin" else "libmqc.so"
    here = os.path.dirname(os.path.abspath(__file__))
    # An in-tree build, which is how this is used before it is installed.
    for up in ("..", "../..", "../../.."):
        yield os.path.join(here, up, "build", name)
    yield name  # let the loader search


def _load():
    tried = []
    for path in _candidates():
        try:
            return ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
        except OSError as exc:
            tried.append(f"  {path}\n    {exc}")
    raise ImportError(
        "could not load libmqc. Build it with `cmake --build build --target "
        f"mqc_shared`, or point {_ENV} at it.\nTried:\n" + "\n".join(tried)
    )


lib = _load()


def _declare(name, restype, argtypes):
    fn = getattr(lib, name)
    fn.restype = restype
    fn.argtypes = argtypes
    return fn


# -- session ---------------------------------------------------------------
session_begin = _declare("mqc_session_begin", _c_int, [])
session_end = _declare("mqc_session_end", _c_int, [])
session_rank = _declare("mqc_session_rank", _c_int, [])
session_size = _declare("mqc_session_size", _c_int, [])
session_is_root = _declare("mqc_session_is_root", _c_int, [])
session_active = _declare("mqc_session_active", _c_int, [])
session_last_error = _declare("mqc_session_last_error", None, [_c_int, _c_str])

# -- system ----------------------------------------------------------------
system_new = _declare("mqc_system_new", _c_ptr, [])
system_free = _declare("mqc_system_free", None, [_c_ptr])
system_set_geometry = _declare(
    "mqc_system_set_geometry",
    _c_int,
    [_c_ptr, _c_int, ctypes.POINTER(_c_int), ctypes.POINTER(_c_double), _c_int, _c_int],
)
system_set_monomers = _declare(
    "mqc_system_set_monomers",
    _c_int,
    [_c_ptr, _c_int, _c_int] + [ctypes.POINTER(_c_int)] * 4,
)
system_set_bonds = _declare(
    "mqc_system_set_bonds", _c_int, [_c_ptr, _c_int] + [ctypes.POINTER(_c_int)] * 4
)
system_n_atoms = _declare("mqc_system_n_atoms", _c_int, [_c_ptr])
system_n_monomers = _declare("mqc_system_n_monomers", _c_int, [_c_ptr])
system_n_bonds = _declare("mqc_system_n_bonds", _c_int, [_c_ptr])
system_bonds_declared = _declare("mqc_system_bonds_declared", _c_int, [_c_ptr])
system_perceive_bonds = _declare("mqc_system_perceive_bonds", _c_int, [_c_ptr, _c_double])
system_count_missing_bonds = _declare(
    "mqc_system_count_missing_bonds", _c_int, [_c_ptr, _c_double]
)
system_auto_monomers = _declare("mqc_system_auto_monomers", _c_int, [_c_ptr, _c_double])
system_last_error = _declare("mqc_system_last_error", None, [_c_int, _c_str])

# -- fragment list ---------------------------------------------------------
fraglist_new = _declare("mqc_fraglist_new", _c_ptr, [])
fraglist_free = _declare("mqc_fraglist_free", None, [_c_ptr])
fraglist_generate = _declare("mqc_fraglist_generate", _c_int, [_c_ptr, _c_int, _c_int])
fraglist_count = _declare("mqc_fraglist_count", _c_long, [_c_ptr])
fraglist_max_level = _declare("mqc_fraglist_max_level", _c_int, [_c_ptr])
fraglist_get = _declare(
    "mqc_fraglist_get", _c_int, [_c_ptr, ctypes.POINTER(_c_int), _c_long, _c_int]
)
fraglist_set = _declare(
    "mqc_fraglist_set", _c_int, [_c_ptr, ctypes.POINTER(_c_int), _c_long, _c_int]
)
fraglist_close_subsets = _declare("mqc_fraglist_close_subsets", _c_int, [_c_ptr])
fraglist_last_error = _declare("mqc_fraglist_last_error", None, [_c_str, _c_int])

# -- run -------------------------------------------------------------------
run = _declare(
    "mqc_run",
    _c_int,
    [_c_ptr, _c_ptr, _c_int, _c_str, _c_int, _c_str, _c_int, ctypes.POINTER(_c_double)],
)
run_last_error = _declare("mqc_run_last_error", None, [_c_int, _c_str])


def message(reader, reversed_args=False):
    """Read a `*_last_error` message back as a str."""
    buf = ctypes.create_string_buffer(512)
    if reversed_args:
        reader(buf, 512)  # fraglist takes (buffer, len)
    else:
        reader(512, buf)
    return buf.value.decode("utf-8", "replace")


def ints(values):
    """A C array of ints from any iterable."""
    values = list(values)
    return (_c_int * max(len(values), 1))(*values)


def doubles(values):
    values = list(values)
    return (_c_double * max(len(values), 1))(*values)
