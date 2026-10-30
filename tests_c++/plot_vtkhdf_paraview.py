"""Read the files of tests_c++/plot_vtkhdf through ParaView's VTKHDF reader (run with pvpython).

    pvpython tests_c++/plot_vtkhdf_paraview.py output/plot_vtkhdf-1-3.vtkhdf [...]

The C++ test checks the files against HDF5; this checks that ParaView itself reads them as a
time series: three time steps, a mesh with points and cells, and a point array "values" whose
last component (the solution u = t x on the unit cube) ranges over [0, t] at time t.
"""

import sys

from paraview.simple import VTKHDFReader

failed = 0
for path in sys.argv[1:]:
    reader = VTKHDFReader(FileName=[path])
    reader.UpdatePipelineInformation()
    times = list(reader.TimestepValues)
    ok = times == [0., 0.5, 1.]
    for t in times:
        reader.UpdatePipeline(t)
        info = reader.GetDataInformation()
        values = reader.PointData["values"]
        n_comp = values.GetNumberOfComponents()
        lo, hi = values.GetRange(n_comp - 1)
        ok &= info.GetNumberOfPoints() > 0 and info.GetNumberOfCells() > 0
        ok &= abs(lo - 0.) < 1e-5 and abs(hi - t) < 1e-5
    print(f"  {path}: {len(times)} steps, {info.GetNumberOfPoints()} points, "
          f"{info.GetNumberOfCells()} cells: {'ok' if ok else 'FAILED'}")
    failed += not ok
sys.exit(1 if failed else 0)
