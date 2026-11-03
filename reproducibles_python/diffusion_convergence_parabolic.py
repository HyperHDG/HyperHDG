#!/usr/bin/env python3
"""Convergence of the parabolic LDG-H diffusion solver (DiffusionParab, Parabolic loop, implicit
Euler, 10^4 time steps to t = 1) on the unit cube in one and two dimensions, degrees 1 to 3, six
refinements: the errors at t = 1 go to output/diffusion_convergence_parabolic.yaml, the
final state of every run is plotted. The classes DiffusionParab_D<dim>_P<degree> come from
modules/diffusion_parabolic.cxx.

    python3 reproducibles_python/diffusion_convergence_parabolic.py [--quick]
"""
import common

NAME = "diffusion_convergence_parabolic"
hdg_module = common.module("diffusion_parabolic")


def diffusion_test(poly_degree, dimension, iteration):
  theta, time_steps = 1., (10 ** 4 if not common.QUICK else 100)
  delta_time = 1 / time_steps
  hdg = hdg_module[f"DiffusionParab_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension,
                                                                  [1., theta, delta_time])
  solution = common.time_stepping(hdg, time_steps, delta_time)
  error = hdg.errors(list(solution), 1.)[0]
  common.report(NAME, error, theta=theta, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"diff_conv_parab{dimension}-{iteration}", time=1.)


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))
