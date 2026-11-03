#!/usr/bin/env python3
"""The parabolic advection solver (AdvectionParab, Parabolic loop, theta scheme with theta = 1/2
and 1) on the unit square with LeVeque's rotating velocity field, 10^3 time steps for one
revolution (t = 2 pi), degrees 0 to 2, 32 to 128 cells per axis: the errors after the revolution
go to output/advection_convergence_rotation.yaml, the final state of every run is
plotted. The classes AdvectionParab_D2_P<degree> come from modules/advection_rotation_parabolic.cxx.

    python3 reproducibles_python/advection_convergence_rotation.py [--quick]
"""
import numpy as np

import common

NAME = "advection_convergence_rotation"
hdg_module = common.module("advection_rotation_parabolic")


def diffusion_test(theta, poly_degree, dimension, iteration):
  time_steps, time_end = (10 ** 3 if not common.QUICK else 10), 2. * np.pi
  delta_time = time_end / time_steps
  hdg = hdg_module[f"AdvectionParab_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension,
                                                                  [0., theta, delta_time])
  solution = common.time_stepping(hdg, time_steps, delta_time, solvers=("gmres", "bicgstab"))
  error = hdg.errors(list(solution), time_end)[0]
  common.report(NAME, error, theta=theta, poly_degree=poly_degree, iteration=iteration)
  common.plot(hdg, solution, f"adv_conv_rotat{theta}-{poly_degree}-{iteration}", time=time_end)


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, theta=[0.5, 1.], poly_degree=range(3), dimension=[2], iteration=range(5, 8) if not common.QUICK else [5])
