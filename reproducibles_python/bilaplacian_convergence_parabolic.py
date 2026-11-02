#!/usr/bin/env python3
"""Convergence of the parabolic LDG-H bilaplacian solver (BilaplacianParab, Parabolic loop,
implicit Euler, 10^4 time steps to t = 1) on the unit cube in one and two dimensions, degrees 1 to
3, six refinements: the errors at t = 1 go to output/bilaplacian_convergence_parabolic.yaml, the
final state of every run is plotted. The classes BilaplacianParab_D<dim>_P<degree> come from
modules/bilaplacian_parabolic.cxx.

    python3 reproducibles_python/bilaplacian_convergence_parabolic.py [--quick]
"""
import common

NAME = "bilaplacian_convergence_parabolic"
hdg_module = common.module("bilaplacian_parabolic")


def bilaplacian_test(poly_degree, dimension, iteration):
  time_steps = 10 ** 4 if not common.QUICK else 100
  delta_time = 1 / time_steps
  hdg = hdg_module[f"BilaplacianParab_D{dimension}_P{poly_degree}"]([2 ** iteration] * dimension,
                                                                    [1., 1., delta_time])
  # as the cython-era script: the data of the last step is not set before the error is taken
  solution = common.time_stepping(hdg, time_steps, delta_time, set_final_data=False)
  error = hdg.errors(list(solution), 1.)[0]
  common.report(NAME, error, poly_degree=poly_degree, dimension=dimension, iteration=iteration)
  common.plot(hdg, solution, f"bil_conv_parab-{dimension}-{iteration}", time=1.)


if __name__ == "__main__":
  common.main(bilaplacian_test, rtol=1e-5, poly_degree=range(1, 4), dimension=range(1, 3), iteration=range(6))
