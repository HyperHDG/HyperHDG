import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_lin_alg
from joblib import Parallel, delayed

import os, sys

try:
	import HyperHDG
except(ImportError, ModuleNotFoundError) as error:
	sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../import")
	import HyperHDG

pol_g = 1

n_shifts = 16
gen_vec = HyperHDG.qmc_methods.generating_vector("lattice-32001-1024-1048576.3600", 100)


zfgen = np.random.default_rng(13)

const = HyperHDG.config()
const.global_loop = "Elliptic"
const.topology = "Cubic<2, 2>"
const.geometry = "UnitCube<2, 2, double>"
const.node_descriptor = "Cubic<2, 2>"
const.local_solver = "Diffusion<2, " + str(pol_g) + "," + str(2*pol_g) + ", TestParametersDiffusionAffine, double, std::vector<double> >"
const.param_t = "std::vector<double>"
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]", "double", "double", "vector[double]"]
const.include_files   = ["reproducibles_python/parameters/diffusion_qmc.hxx"]
const.debug_mode      = False

PyDP = HyperHDG.include(const)
HDG_wrapper = PyDP([100, 100])

system_size = HDG_wrapper.size_of_system()

with open("output/kvgr.txt", "w") as f:
	f.write("n_qmc_points\terr_u\terr_q\n")

def solve(k, gen_vec, shift, n_qmc_points, HDG_wrapper):
	print("Quadrature point " + str(k+1) + " of " + str(n_qmc_points))
	system_size = HDG_wrapper.size_of_system()
	arr = HyperHDG.qmc_methods.get_quadrature_point(gen_vec, shift, k+1, n_qmc_points) - 0.5
	vectorRHS = np.multiply( HDG_wrapper.residual_flux(np.zeros(system_size), arr), -1. )
	col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(arr)
	A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))
	[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, rtol=1e-11)
	return HDG_wrapper.mean(vectorSolution, arr)

for i in range(2, 13):
	n_qmc_points = 2**i
	print("Quadrature points: " + str(n_qmc_points))
	mean_sa = 0
	expec = np.zeros((n_shifts, 2))
	for m in range(n_shifts):
		print("Shift " + str(m+1) + " of " + str(n_shifts))
		shift = zfgen.random(100)
		means = np.array(Parallel(n_jobs=-2, prefer="threads")( delayed(solve) (k, gen_vec, shift, n_qmc_points, HDG_wrapper) for k in range(n_qmc_points)))
		print(means)
		expec[m, 0] = np.mean(means[:, -1])
		expec[m, 1] = np.mean(np.linalg.norm(means[:, :-1]))
	print(expec)
	mean_sa = np.mean(expec, axis=0)
	print(mean_sa)
	print(expec - mean_sa)
	err = np.linalg.norm(mean_sa - expec, axis=0) / np.sqrt(n_shifts - 1) / np.sqrt(n_shifts)
	print(err)
	with open("output/kvgr.txt", "a") as f:
		f.write(str(n_qmc_points) + "\t" + str(err[0]) + "\t" + str(err[1]) + "\n")

