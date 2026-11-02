#!/usr/bin/env python3
"""Convergence of the parabolic advection solver (AdvectionParab, Parabolic loop, theta scheme with
theta = 1/2 and 1, 10^4 time steps to t = 1) on the unit cube in one and two dimensions, degrees 1
to 3, six refinements: the errors at t = 1 go to
output/advection_convergence_parabolic.yaml, the final state of every run is plotted.
The classes AdvectionParab_D<dim>_P<degree> come from modules/advection_parabolic.cxx.

    python3 reproducibles_python/advection_convergence_parabolic.py [--quick]
"""
import common

NAME = "advection_convergence_parabolic"
hdg_module = common.module("advection_parabolic")


def diffusion_test(theta, poly_degree, dimension, iteration):
  time_steps, time_end = (10 ** 4 if not common.QUICK else 100), 1.
  delta_time = time_end / time_steps
  hdg = hdg_module[f"AdvectionParab_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension,
                                                                  [0., theta, delta_time])
  solution = common.time_stepping(hdg, time_steps, delta_time, solvers=("gmres", "bicgstab"))
  error = hdg.errors(list(solution), time_end)[0]
  common.report(NAME, error, theta=theta, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"adv_conv_parab{dimension}-{iteration}", time=time_end)


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, theta=[0.5, 1.], poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))
