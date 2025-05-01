import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_lin_alg
from scipy.sparse.linalg import LinearOperator

import os, sys

try:
	import HyperHDG
except(ImportError, ModuleNotFoundError) as error:
	sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
	import HyperHDG

pol_g = 1

n_shifts = 16
gen_vec = np.loadtxt("reproducibles_python/parameters/lattice-32001-1024-1048576.3600", usecols=1, max_rows=100)

def get_qp(gen, shift, k, n_qmc_points):
	asg = k * gen / n_qmc_points
	asg = asg - np.floor(asg) + shift
	asg = asg - np.floor(asg) - 0.5
	return asg

zfgen = np.random.default_rng(27)

const = HyperHDG.config()
const.global_loop = "Elliptic"
const.topology = "Cubic<2, 2>"
const.geometry = "UnitCube<2, 2, double>"
const.node_descriptor = "Cubic<2, 2>"
const.local_solver = "Diffusion<2, " + str(pol_g) + "," + str(2*pol_g) + ", TestParametersDiffusionAffine, double, std::vector<double> >"
const.param_t = "std::vector<double>"
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]"]
const.include_files   = ["reproducibles_python/parameters/diffusion_qmc.hxx"]
const.debug_mode      = False

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP([100, 100])

system_size = HDG_wrapper.size_of_system()
mean = np.zeros(system_size)
mean_sa = np.zeros(system_size)

with open("output/kvgr.txt", "w") as f:
	f.write("n_qmc_points\tsq_err\n")
	f.close()

for i in range(2, 6):
	n_qmc_points = 2**i
	print("Quadraturpunkte: " + str(n_qmc_points))
	mean_sa = np.zeros(system_size)
	for m in range(n_shifts):
		print("Shift " + str(m+1) + " von " + str(n_shifts))
		shift = zfgen.random(100)
		mean = np.zeros(system_size)
		for k in range(n_qmc_points):
			print("Punkt " + str(k+1) + " von " + str(n_qmc_points))
			arr = get_qp(gen_vec, shift, k+1, n_qmc_points)
			vectorRHS = np.multiply( HDG_wrapper.residual_flux(np.zeros(system_size), arr), -1. )
			col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(arr)
			A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))
			[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, tol=1e-11)
			mean += vectorSolution / n_qmc_points
		mean_sa += mean / n_shifts
	with open("output/kvgr.txt", "a") as f:
		f.write(str(n_qmc_points) + "\t" + str(np.linalg.norm(mean_sa - mean)) + "\n")
		f.close()

