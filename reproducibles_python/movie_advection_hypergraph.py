#!/usr/bin/env python3
"""A movie of the parabolic advection solver (AdvectionParab, Parabolic loop, theta scheme with
theta = 1/2 and 1, 10^3 time steps to t = 5, GMRES tolerance 1e-8) on the hypergraph of
two-dimensional hyperedges in domains/leVeque_hg.geo, every hyperedge split 32 times, degrees 1 to
3: the state is plotted at the start and after every step, the errors at t = 5 go to
output/movie_advection_hypergraph.yaml.
The classes AdvectionParab_D2_P<degree> come from modules/advection_movie_file_parabolic.cxx.

    python3 reproducibles_python/movie_advection_hypergraph.py [--quick]
"""
import common

NAME = "movie_advection_hypergraph"
hdg_module = common.module("advection_movie_file_parabolic")
DOMAIN = str(common.ROOT / "domains" / "leVeque_hg.geo")


def diffusion_test(theta, poly_degree, refinement):
  time_steps, time_end = (10 ** 3 if not common.QUICK else 20), 5
  delta_time = time_end / time_steps
  hdg = hdg_module[f"AdvectionParab_D2_P{poly_degree}"](DOMAIN, [0., theta, delta_time])
  hdg.set_refinement(2 ** refinement)
  filename = f"leVeque_hyg{theta}-{poly_degree}-{refinement}"
  solution = common.time_stepping(
    hdg, time_steps, delta_time, tol=1e-8, solvers=("gmres", "bicgstab"),
    on_start=lambda x: common.plot(hdg, x, filename, time=time_end),
    on_step=lambda step, x: hdg.plot_solution(list(x), time_end))
  error = hdg.errors(list(solution), time_end)[0]
  common.report(NAME, error, theta=theta, poly_degree=poly_degree, refinement=refinement)


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-4, theta=[0.5, 1.], poly_degree=range(1, 4), refinement=range(5, 6))
