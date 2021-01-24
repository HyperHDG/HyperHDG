from __future__ import print_function
import numpy, os, sys
import scipy.sparse.linalg as sp_lin_alg

try:
  import HyperHDG
except (ImportError, ModuleNotFoundError) as error:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
  import HyperHDG

# Importing done

poly_degree = 1
hyEdge_dim  = 1
cube_dim    = 2
refinement  = 1
debug_mode  = True
  
hdg_config                 = HyperHDG.config()
hdg_config.global_loop     = "Elliptic"
hdg_config.topology        = f'Cubic<{hyEdge_dim},{cube_dim}>'
hdg_config.geometry        = f'UnitCube<{hyEdge_dim},{cube_dim},double>'
hdg_config.node_descriptor = f'Cubic<{hyEdge_dim},{cube_dim}>'
hdg_config.local_solver    = "Diffusion<" + str(hyEdge_dim) + "," + str(poly_degree) + "," + str(
  2*poly_degree) + ",HG<" + str(hyEdge_dim) + ">::DiffusionElliptic,double>"
hdg_config.include_files   = ["examples/parameters/diffusion.hxx"]
hdg_config.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
hdg_config.debug_mode      = debug_mode

hyperHDGClass = HyperHDG.include(hdg_config)
HDGObject   = hyperHDGClass ( [2 ** refinement] * cube_dim )

vectorRHS = numpy.multiply( HDGObject.residual_flux(HDGObject.zero_vector()), -1. )

system_size = HDGObject.size_of_system()
A = sp_lin_alg.LinearOperator((system_size,system_size), matvec=HDGObject.trace_to_flux)

[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-13)
if num_iter != 0:
  print("CG solver failed with a total number of ", num_iter, "iterations.")
  raise RuntimeError("Linear solvers did not converge!")

print("Error: ", HDGObject.errors(vectorSolution)[0])

HDGObject.plot_option( "fileName" , "diffusion_elliptic_py" )
HDGObject.plot_option( "printFileNumber" , "false" )
HDGObject.plot_option( "scale" , "0.95" )
HDGObject.plot_solution(vectorSolution)
