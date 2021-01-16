from __future__ import print_function
import numpy, os, sys
import scipy.sparse.linalg as sp_lin_alg

poly_degree = 1
hyEdge_dim  = 1
space_dim   = 2
refinement  = 1
debug_mode  = True
  
try:
  import HyperHDG
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
  import HyperHDG
  
const                 = HyperHDG.config()
const.global_loop     = "Elliptic"
const.topology        = "Cubic<" + str(hyEdge_dim) + "," + str(space_dim) + ">"
const.geometry        = "UnitCube<" + str(hyEdge_dim) + "," + str(space_dim) + ",double>"
const.node_descriptor = "Cubic<" + str(hyEdge_dim) + "," + str(space_dim) + ">"
const.local_solver    = "Diffusion<" + str(hyEdge_dim) + "," + str(poly_degree) + "," + str(
  2*poly_degree) + ",HG<" + str(hyEdge_dim) + ">::DiffusionElliptic,double>"
const.include_files   = ["examples/parameters/diffusion.hxx"]
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
const.debug_mode      = debug_mode

hyperHDG    = HyperHDG.include(const)
HDG_wrapper = hyperHDG( [2 ** refinement] * space_dim )

vectorRHS = numpy.multiply( HDG_wrapper.residual_flux(HDG_wrapper.zero_vector()), -1. )

system_size = HDG_wrapper.size_of_system()
A = sp_lin_alg.LinearOperator((system_size,system_size), matvec=HDG_wrapper.trace_to_flux)

[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)
if num_iter != 0:
  print("CG solver failed with a total number of ", num_iter, "iterations.")
  raise RuntimeError("Linear solvers did not converge!")

print("Error: ", HDG_wrapper.errors(vectorSolution)[0])

HDG_wrapper.plot_option( "fileName" , "diffusion_elliptic_py" )
HDG_wrapper.plot_option( "printFileNumber" , "false" )
HDG_wrapper.plot_option( "scale" , "0.95" )
HDG_wrapper.plot_solution(vectorSolution)
