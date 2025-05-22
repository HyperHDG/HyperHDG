import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sp_lin_alg
from scipy.stats import norm as sp_norm
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

zfgen = np.random.default_rng(27)

def phi_affine(arr):
	return arr - 0.5

def phi_lognormal(arr):
	return sp_norm.ppf(arr)

const = HyperHDG.config()
const.global_loop = "Elliptic"
const.topology = "Cubic<2, 2>"
const.geometry = "UnitCube<2, 2, double>"
const.node_descriptor = "Cubic<2, 2>"
const.param_t = "std::vector<double>"
const.cython_replacements = ["vector[unsigned int]", "vector[unsigned int]", "double", "double", "vector[double]"]
const.include_files   = ["reproducibles_python/parameters/diffusion_qmc.hxx"]
const.debug_mode      = False


def solve(k, gen_vec, shift, n_qmc_points, morph):
	HDG_wrapper = PyDP([100, 100])

	system_size = HDG_wrapper.size_of_system()
	
	print("Quadrature point " + str(k+1) + " of " + str(n_qmc_points))
	system_size = HDG_wrapper.size_of_system()
	arr = morph(HyperHDG.qmc_methods.get_quadrature_point(gen_vec, shift, k+1, n_qmc_points))
	vectorRHS = np.multiply( HDG_wrapper.residual_flux(np.zeros(system_size), arr), -1. )
	col_ind, row_ind, vals = HDG_wrapper.sparse_stiff_mat(arr)
	A = sp.csr_matrix((vals, (row_ind,col_ind)), shape=(system_size,system_size))
	[vectorSolution, num_iter] = sp_lin_alg.cg(A, vectorRHS, rtol=1e-11)
	return [HDG_wrapper.mean(vectorSolution, arr)]

def get_err(gen_vec, n_qmc_points, morph):
	mean_sa = 0
	expec = np.zeros((n_shifts, 4))
	for m in range(n_shifts):
		print("Shift " + str(m+1) + " of " + str(n_shifts))
		shift = zfgen.random(100)
		asg = Parallel(n_jobs=32)( delayed(solve) (k, gen_vec, shift, n_qmc_points, phi_affine) for k in range(n_qmc_points))
		means = np.array([run[0] for run in asg])
		expec[m, 0] = np.mean(np.sqrt(means[:, 0]))
		expec[m, 1] = np.mean(np.sqrt(means[:, 1]))
		expec[m, 2] = np.mean(np.sqrt(means[:, 2]))
		#scelet_v = np.array([run[1] for run in asg])
		#expec[m, 3] = np.mean(np.linalg.norm(scelet_v, axis=1))
		print(expec[m])
	mean_sa = np.mean(expec, axis=0)
	err = np.linalg.norm(mean_sa - expec, axis=0) / np.sqrt(n_shifts - 1) / np.sqrt(n_shifts)
	print(err)
	return err


#affine case

const.local_solver = "Diffusion<2, " + str(pol_g) + "," + str(2*pol_g) + ", TestParametersDiffusionAffine, double, std::vector<double> >"
PyDP = HyperHDG.include(const)

with open("output/kvgr_affine.txt", "w") as f:
	f.write("n_qmc_points\terr_u\terr_q\terr_g\terr_s\n")

for i in range(2, 12):
	n_qmc_points = 2**i
	err = get_err(gen_vec, n_qmc_points, phi_affine)
	with open("output/kvgr_affine.txt", "a") as f:
		f.write(str(n_qmc_points) + "\t" + str(err[0]) + "\t" + str(err[1]) + "\t" + str(err[2]) + "\t" + str(err[3]) + "\n")


#lognormal case

const.local_solver = "Diffusion<2, " + str(pol_g) + "," + str(2*pol_g) + ", TestParametersDiffusionLognormal, double, std::vector<double> >"
PyDP = HyperHDG.include(const)

with open("output/kvgr_lognormal.txt", "w") as f:
	f.write("n_qmc_points\terr_u\terr_q\terr_g\terr_s\n")

for i in range(2, 12):
	n_qmc_points = 2**i
	err = get_err(gen_vec, n_qmc_points, phi_lognormal)
	with open("output/kvgr_lognormal.txt", "a") as f:
		f.write(str(n_qmc_points) + "\t" + str(err[0]) + "\t" + str(err[1]) + "\t" + str(err[2]) + "\t" + str(err[3]) + "\n")


#Gevrey affine case

const.local_solver = "Diffusion<2, " + str(pol_g) + "," + str(2*pol_g) + ", TestParametersDiffusionGevreyAffine, double, std::vector<double> >"
PyDP = HyperHDG.include(const)

with open("output/kvgr_gevrey_affine.txt", "w") as f:
	f.write("n_qmc_points\terr_u\terr_q\terr_g\terr_s\n")

for i in range(2, 12):
	n_qmc_points = 2**i
	err = get_err(gen_vec, n_qmc_points, phi_affine)
	with open("output/kvgr_gevrey_affine.txt", "a") as f:
		f.write(str(n_qmc_points) + "\t" + str(err[0]) + "\t" + str(err[1]) + "\t" + str(err[2]) + "\t" + str(err[3]) + "\n")


#Gevrey lognormal case

const.local_solver = "Diffusion<2, " + str(pol_g) + "," + str(2*pol_g) + ", TestParametersDiffusionGevreyLognormal, double, std::vector<double> >"
PyDP = HyperHDG.include(const)

with open("output/kvgr_gevrey_lognormal.txt", "w") as f:
	f.write("n_qmc_points\terr_u\terr_q\terr_g\terr_s\n")

for i in range(2, 12):
	n_qmc_points = 2**i
	err = get_err(gen_vec, n_qmc_points, phi_lognormal)
	with open("output/kvgr_gevrey_lognormal.txt", "a") as f:
		f.write(str(n_qmc_points) + "\t" + str(err[0]) + "\t" + str(err[1]) + "\t" + str(err[2]) + "\t" + str(err[3]) + "\n")
