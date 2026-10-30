# Command-line tools

Built with the library (`cmake --build build`), executables land in `build/tools/`.

## geo2h5

Converts a `.geo` domain file into its HDF5 counterpart `.geo.h5`
(format: `include/HyperHDG/read_domain_hdf5.hxx`). Needs HyperHDG configured with HDF5.

```sh
build/tools/geo2h5 domains/cross.geo                 # -> domains/cross.geo.h5
build/tools/geo2h5 domains/cross.geo out/cross.geo.h5
build/tools/geo2h5 domains/cross.geo --no-vtkhdf     # without the ParaView view
```

The hyperedge and space dimensions are read from the file header; the combinations compiled in
are listed in the `switch` at the end of `geo2h5.cxx`.

## Inspecting `.h5` files

HDF5 ships `h5dump` and `h5ls` (packages `hdf5` or `hdf5-tools`; also in a conda environment with
`hdf5`). The useful calls, none of which print the data unless asked:

```sh
h5dump -n  domains/cross.geo.h5               # list groups, datasets and links
h5dump -H  domains/cross.geo.h5               # shapes and types of everything, no data
h5dump -A  domains/cross.geo.h5               # attributes (space_dim, hyperedge_dim, ...)
h5dump -d /domain/edges domains/cross.geo.h5  # one dataset with its data
h5dump -p -d /domain/points domains/cross.geo.h5   # plus chunking and compression
h5ls -r    domains/cross.geo.h5               # one line per object, with shapes
```

From Python, h5py does the same:

```python
import h5py
with h5py.File("domains/cross.geo.h5") as f:
    f.visititems(lambda name, obj: print(name, getattr(obj, "shape", "")))
    points, edges = f["domain/points"][:], f["domain/edges"][:]
```

ParaView reads the `/VTKHDF` view as an unstructured grid. It picks readers by file extension and
does not know `.geo.h5`, so choose "VTKHDF Reader" when it asks, or symlink the file to a `.vtkhdf`
name.
