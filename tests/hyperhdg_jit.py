#!/usr/bin/env python3
"""python/hyperhdg.py: compiling module sources just in time.

- A module source is compiled for the running interpreter, and the extension module is found and
  imported; its members are reached by attribute and by subscript.
- Variants of a module (hash_name=True) get build directories of their own, so switching back to
  one does not rebuild it, and load imports an edited variant within the same session.

The modules are trivial (no HyperHDG headers), so with ccache they compile from the cache even
when HyperHDG changes; that a HyperHDG module builds just in time is shown by
examples/diffusion_elliptic_nb_jit.py, the NumPy bindings are tested by tests/bindings.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))
import hyperhdg  # noqa: E402

failures = 0


def check(ok, what):
  global failures
  print(f"{'ok  ' if ok else 'FAIL'}  {what}")
  failures += not ok


def code(value):
  return ("#include <nanobind/nanobind.h>\n"
          "NB_MODULE(hyperhdg_jit_test, m) { m.def(\"value\", []() { return %d; }); }\n" % value)


# ---- compiling, finding and importing
module = hyperhdg.load(code(1))
check(module.value() == 1, "the module is compiled, found and imported")
check(module["value"] is module.value, "members by attribute and by subscript")

# ---- variants: separate build directories, no rebuild when switching back
paths = [Path(hyperhdg.compile(code(value), hash_name=True)) for value in (1, 2)]
mtime = paths[0].stat().st_mtime_ns
again = Path(hyperhdg.compile(code(1), hash_name=True))
check(paths[0] != paths[1] and paths[0].parent != paths[1].parent,
      "two variants are built in separate directories")
check(again == paths[0] and again.stat().st_mtime_ns == mtime,
      "switching back to a variant does not rebuild it")
check(hyperhdg.load(code(2)).value() == 2 and module.value() == 1,
      "load imports an edited variant next to the first one")

sys.exit(1 if failures else 0)
