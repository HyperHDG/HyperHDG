//using namespace HyperHDG;
/**
 * \page exdiff Example Python code for diffusion

Using HyperHDG in a Python script can be done by exploiting the script
`cython_import.py`. This way is exemplary discussed for the an
elliptic diffusion problem with an example application. The script can
be found in HyperHDG's file `examples/diffusion_elliptic.py`.

We start with importing the print function, and the packages numpy,
os, sys, and scipy's linear algebra tools. These packages will be used
in the remainder of the Python script.

\dontinclude diffusion_elliptic.py
\until sp_lin_alg

Then, we import our own importer for C++ code based on cython.

\until done

Next, we define some basic parameters such as the polynomial degree of
both, the (discontinuous) polynomials which live on the skeletal, and
the polynomials which are utilized by the local solvers, the dimension
of a hyperedge, the dimension of the surrounding space, the refinement
level of the domain, and whether we want the C++ code to be compiled
in debug or release mode.

\until debug

Now we configure the discretization we want to work on as a textual
representation, which is then fed into `cython_import`.

\until debug

This text is compiled into C++ code and then to a linkable object file, which creates a Python class. This in turn is used to create an object
\todo Why don't we construct in cython right away?
\todo explain arguments

\until wrapper

Now, we can start dealing with the actual problem. To do so, we create the global right hand side by evaluating the residual of A x - b for x = 0 and multiplying the result by -1.

\until multiply

Python's Scipy package provides solvers for matrix-free situations. These solver require the definition of a LinearOperator, which is done in two lines:

\until sp_lin_alg

Solving the linear system of equations can then be done using the CG method in one line. The remaining lines check that the CG method has converged.

\until RuntimeError

In the case that the analytical solution is known, we can evaluate the error and print it using:

\until Error

Additionally, we can plot the calculated solution after setting some options found in PlotOptions:

\until plot_solution

 */
