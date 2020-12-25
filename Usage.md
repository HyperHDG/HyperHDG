# How to write your own program using HyperHDG?


## Write a Python script

Using HyperHDG in a Python script can be done by exploiting the script `cython_import`. This way is
exemplary discussed for the an elliptic diffusion problem with an example application. The script
can be found in HyperHDG's file `examples/diffusion_elliptic.py`.

We start with importing the print function, and the packages `numpy`, `os`, `sys`, and `scipy`'s
linear algebra tools. These packages will be used in the remainder of the Python script.

```
from __future__ import print_function
import numpy, os, sys
import scipy.sparse.linalg as sp_lin_alg
```

Next, we define some basic parameters such as the polynomial degree of both, the (discontinuous)
polynomials which live on the skeletal, and the polynomials which are utilized by the local solvers,
the dimension of a hyperedge `hyEdge_dim`, the dimension of the surrounding space `space_dim`, the
refinement level of the domain, and whether we want the C++ code to be compiled in debug or release
mode.

```
poly_degree = 1
hyEdge_dim  = 1
space_dim   = 2
refinement  = 1
debug_mode  = True
```

We have not used HyperHDG so far. This will be changed in the next lines which try to import the
interface `cython_import`. If this is not successful, we add the path to `cython_import`'s
directory, which is the main directory of the repository, manually. Thus, the expression 
`os.path.dirname( os.path.abspath( __file__ ) ) + "/.."` might be changed to the path of
HyperHDG in self-written scripts.

```
try:
  import cython_import
except ImportError as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
  import cython_import
```

Having imported HyperHDG, we want to configure the possible problem, HyperHDG should deal with. That
is why, we create a `hyperhdg_constructor` and set its properties to define our problem. Our problem
is a standard elliptic problem, for which an elliptic global loop can be used. To inform ourselves
about the features of an elliptic loop and those of possible other loops, the [Doxygen](
https://andreasrupp.github.io/HyperHDG_pages/doxygen) pages may be used. Our problem will be defined
on a `space_dim` dimensional (refined) unit cube's `hyEdge_dim` dimensional elements. Thus, we need
a cubic topology, a cubic node--descriptor, and the geometry of a unit cube. Moreover, the local
solver to be chosen is our diffusion solver. Again, the meaning of the respective template arguments
can be found in the [Doxygen](https://andreasrupp.github.io/HyperHDG_pages/doxygen). Here, 
`poly_degree` refers to the order of the polynomial space, `2 * poly_degree` refers to the order of
the quadrature rule, the next parameter refers to a `struct` in which the utilized parameters (such
as right-hand side, diffusivity, ...) are defined, and `double` is the used floating point 
arithmetic of the local solvers. The file in which the parameters `struct` is defined needs to be
known to HyperHDG. Thus, it is added to a list of (additional) include files. The other files are
automatically found by `cython_import`. The cython replacements define the constructor arguments
(also denoted `constructor_values` in the [Doxygen](
https://andreasrupp.github.io/HyperHDG_pages/doxygen)) of the topology and the geometry, and the
debug mode is also transmitted to `cython_import`.

```
const                 = cython_import.hyperhdg_constructor()
const.global_loop     = "Elliptic"
const.topology        = "Cubic<" + str(hyEdge_dim) + "," + str(space_dim) + ">"
const.geometry        = "UnitCube<" + str(hyEdge_dim) + "," + str(space_dim) + ",double>"
const.node_descriptor = "Cubic<" + str(hyEdge_dim) + "," + str(space_dim) + ">"
const.local_solver    = "Diffusion<" + str(hyEdge_dim) + "," + str(poly_degree) + "," + str(
  2*poly_degree) + ",HG<" + str(hyEdge_dim) + ">::DiffusionElliptic,double>"
const.include_files   = ["examples/parameters/diffusion.hxx"]
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
const.debug_mode      = debug_mode
```

With the statical configuration (the one that is needed for compilation of the C++ files) done, we
can now define the class `hyperHDG` and create an element of this class by using its constructor.
Note that it automatically fills the constructors for geometry and local solver with appropriate
values. Nonetheless, the constructor might also comprise up to three elements. In this case, two
vectors and an positive number (which represents the penalty in the LDG-H method). 

```
hyperHDG    = cython_import.cython_import(const)
HDG_wrapper = hyperHDG( [2 ** refinement] * space_dim )
```

Now, we can start dealing with the actual problem. To do so, we create the global right hand side by
evaluating the residual of `A x - b` for `x = 0` and multiplying the result by `-1`.

```
vectorRHS = numpy.multiply( HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector()), -1. )
```

The matrix `A` will not be calculated by HyperHDG. Instead, the operator mapping `x` to `A x` is
realized in a matrix-free fashion. Python's `scipy` package provides solvers for matrix-free
situations. These solver require the definition of a `LinearOperator`, which is done in two lines:

```
system_size = HDG_wrapper.size_of_system()
A = sp_lin_alg.LinearOperator((system_size,system_size), matvec=HDG_wrapper.matrix_vector_multiply)
```

Solving the linear system of equations can then be done using the CG method in one line. The
remaining lines check that the CG method has converged.

```
[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)
if num_iter != 0:
  print("CG solver failed with a total number of ", num_iter, "iterations.")
  raise RuntimeError("Linear solvers did not converge!")
```

In the case that the analytical solution is known, we can evaluate the error and print it using:

```
print("Error: ", HDG_wrapper.calculate_L2_error(vectorSolution))
```

Additionally, we can plot the calculated solution after setting some plot options whose meaning can
be found in the [Doxygen](https://andreasrupp.github.io/HyperHDG_pages/doxygen).

```
HDG_wrapper.plot_option( "fileName" , "diffusion_elliptic_py" )
HDG_wrapper.plot_option( "printFileNumber" , "false" )
HDG_wrapper.plot_option( "scale" , "0.95" )
HDG_wrapper.plot_solution(vectorSolution)
```


## Write a C++ program

Alternatively, one can use HyperHDG as C++ library. Please note that for this purpose, the `include`
directory of HyperHDG as well as the `include` directories of the submodules need to be added to the
compiler's include paths. The flag `-DNDEBUG` compiles in release mode, while the code is created in
debug mode if the flag is missing. Moreover, HyperHDG needs to be linked to LAPACK, compiled with a
standard of C++17, and the compiler needs to be able to find the parameters file in this example.
Please check that your compiler meets the standards of HyperHDG. Compilers that are tested against
in pushes to the `master` branch and therefore should be safe can be found in the `Makefile`.

The code of the C++ program that does exactly the same as the aforementioned Python script can be
found in the file `examples/diffusion_elliptic.cxx`.

Since for C++, we do not have `cython_import` that finds all the necessary files, we need to include
them by hand. Since, we want to conduct the same example again, we have to include the geometry and
the topology of a unit cube, as well as the elliptic global loop and the diffusion solver.
Additionally, we need a sparse linear algebra providing the CG solver. Last, but not least the
parameters need to be included.

```
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <examples/parameters/diffusion.hxx>
```

This has prepared us to define the `main` function which will surround the rest of the code

```
int main()
{
```

and starts with defining the very same parameters as before.

```
  const unsigned int poly_degree = 1;
  const unsigned int hyEdge_dim = 1;
  const unsigned int space_dim = 2;
  const unsigned int refinement = 1;
```

Next, we define the problem class, where the different components have exactly the same meaning as
in the `hyperhdg_constructor` of the Python script. The only difference consists in the fact that
this class can immediately be instantiated, which has not been the case in Python, where
`cython_import` had to build it first.

```
  GlobalLoop::Elliptic<Topology::Cubic<hyEdge_dim, space_dim>,
                       Geometry::UnitCube<hyEdge_dim, space_dim>,
                       NodeDescriptor::Cubic<hyEdge_dim, space_dim>,
                       LocalSolver::Diffusion<hyEdge_dim, poly_degree, 2 * poly_degree,
                                              HG<hyEdge_dim>::DiffusionElliptic, double> >
    HDG_wrapper(std::vector<unsigned int>(space_dim, 1 << refinement));
```

The remainder of the program also works analogously to the Python script, i.e., we first define the
right-hand side as negative, homogeneous residuum,

```
  std::vector<double> vectorRHS = HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector());
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)
    vectorRHS[i] *= -1.;
```

solve the global system of equations using the CG method,

```
  std::vector<double> vectorSolution = SparseLA::conjugate_gradient(vectorRHS, HDG_wrapper);
```

and print the L^2 error.

```
  std::cout << "Error: " << HDG_wrapper.calculate_L2_error(vectorSolution) << std::endl;
```

At the end, we plot the solution

```
  HDG_wrapper.plot_option("fileName", "diffusion_elliptic_c++");
  HDG_wrapper.plot_option("printFileNumber", "false");
  HDG_wrapper.plot_option("scale", "0.95");
  HDG_wrapper.plot_solution(vectorSolution);
```

and end the program returning 0.

```
  return 0;
}
```
