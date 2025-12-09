#!/usr/bin/env python3

import numpy as np
import jprecond

# '<' for little-endian
# u8 : uint64_t
# S8 : char[8]

DataTableType = np.dtype([
  ('name', 'S8'),
  ('offset', '<u8'),
  ('size', '<u8'),
])

GeoBinHeaderType = np.dtype([
    ('magic', 'S8'),
    ('space_dim', '<u8'),
    ('hyperedge_dim', '<u8'),
    ('n_points', '<u8'),
    ('n_hypernodes', '<u8'),
    ('n_hyperedges', '<u8'),
])

DomainsHeaderType = np.dtype([
  ('magic', 'S8'),
  ('idsize', '<u8'),
  ('n_domains', '<u8'),
  ('tables', DataTableType, (2,))
])

FloatType = np.dtype("float64")
IdType = np.dtype("<u4")


def read_network_points(path):
  with open(path, "rb") as decom_file:
    header_without_tables = np.frombuffer(
        decom_file.read(GeoBinHeaderType.itemsize),
        dtype=GeoBinHeaderType,
        offset=0,
        count=1
    )[0]
    tables = np.frombuffer(
        decom_file.read(5*DataTableType.itemsize),
        dtype=DataTableType,
        count=5
    )
    size = tables[0]["size"]
    network_points = np.frombuffer(
        decom_file.read(size),
        dtype=FloatType,
        count=int(size/FloatType.itemsize),
    ).reshape(-1, header_without_tables["space_dim"])
    return network_points

def read_domains(path):
  with open(path, "rb") as file:
    header = np.frombuffer(
      file.read(DomainsHeaderType.itemsize),
      dtype=DomainsHeaderType,
    )[0]
    ioffsets = np.frombuffer(
      file.read(header["tables"][0]["size"]),
      dtype=IdType,
    )
    all_domains = np.frombuffer(
      file.read(header["tables"][1]["size"]),
      dtype=IdType,
    )
    domains = jprecond.Domains(ioffsets, all_domains)
    return domains
