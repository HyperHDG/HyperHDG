# External documentation {#external}

HyperHDG builds on three libraries, included as git submodules:

- [TPP](tpp/index.html) (`submodules/tensor_product_polynomials.git`): tensor product polynomials,
  i.e. the shape functions and quadrature rules of the local solvers.
- [TPCC](tpcc/index.html) (`submodules/tensor_product_chain_complex.git`): tensor product chain
  complexes, i.e. the cubic topologies and their boundaries.
- [nanobind](https://nanobind.readthedocs.io) (`submodules/nanobind.git`): the Python bindings
  (`include/HyperHDG/bind_python.hxx`, `python/hyperhdg.py`).

The pages of TPP and TPCC are built together with these; types of TPP and TPCC in HyperHDG's pages
link to them.
