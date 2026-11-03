#!/usr/bin/env python3
"""Convert a .geo domain file into its HDF5 counterpart .geo.h5.

    python python/geo2h5.py domains/cross.geo [out.geo.h5] [--no-vtkhdf]

The format is documented in include/HyperHDG/read_domain_hdf5.hxx. The dimensions are template
parameters of read_domain, so they are read from the file header and the matching instantiation
is compiled on demand (once per combination, cached afterwards) through hyperhdg.load.
"""

import argparse
import os
import re

import hyperhdg  # python/hyperhdg.py, next to this script

MODULE = """\
#include <HyperHDG/read_domain.hxx>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

NB_MODULE(geo2h5, m) {{
  m.def("convert", [](const std::string& input, const std::string& output, bool vtkhdf) {{
    const auto domain = read_domain<{hyEdge_dim}, {space_dim}>(input);
    write_domain_hdf5(output, domain, vtkhdf);
    return std::make_tuple(domain.n_points, domain.n_hyNodes, domain.n_hyEdges,
                           domain.n_properties);
  }});
}}
"""


def header_value(path, keyword):
    """Value of a 'Keyword = N;' line of a .geo file."""
    with open(path) as file:
        return int(re.search(keyword + r"\s*=\s*(\d+)", file.read()).group(1))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("input", help=".geo file")
    parser.add_argument("output", nargs="?", help=".geo.h5 file (default: input + .h5)")
    parser.add_argument("--no-vtkhdf", action="store_true", help="omit the ParaView view")
    args = parser.parse_args()
    output = args.output or args.input + ".h5"

    hyEdge_dim = header_value(args.input, "HyperEdge_Dim")
    space_dim = header_value(args.input, "Space_Dim")
    module = hyperhdg.load(MODULE.format(hyEdge_dim=hyEdge_dim, space_dim=space_dim))
    n_points, n_hyNodes, n_hyEdges, n_properties = module.convert(
        args.input, output, not args.no_vtkhdf)
    print(f"{args.input} -> {output}: {n_points} points, {n_hyNodes} hypernodes, {n_hyEdges} "
          f"hyperedges of dimension {hyEdge_dim} in {space_dim}d, {n_properties} properties")


if __name__ == "__main__":
    main()
