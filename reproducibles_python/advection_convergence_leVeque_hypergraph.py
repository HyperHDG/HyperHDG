#!/usr/bin/env python3
"""The parabolic advection solver (AdvectionParab, Parabolic loop, theta scheme with theta = 1/2 and
1, 10^3 time steps to t = 5) with LeVeque's rotating velocity field on the hypergraph of
two-dimensional hyperedges in domains/leVeque_hg.geo, every hyperedge split 2^refinement times,
degrees 0 to 2, refinements 5 to 7: the errors at t = 5 go to
output/advection_convergence_leVeque_hypergraph.yaml, the state is plotted at the start and every 20 steps. The
classes AdvectionParab_D2_P<degree> come from modules/advection_leveque_file_parabolic.cxx.

    python3 reproducibles_python/advection_convergence_leVeque_hypergraph.py [--quick]
"""
import common

NAME = "advection_convergence_leVeque_hypergraph"
hdg_module = common.module("advection_leveque_file_parabolic")
DOMAIN = str(common.ROOT / "domains" / "leVeque_hg.geo")


def diffusion_test(theta, poly_degree, refinement):
  time_steps, time_end, output_interval = (10 ** 3 if not common.QUICK else 20), 5, 20
  delta_time = time_end / time_steps
  hdg = hdg_module[f"AdvectionParab_D2_P{poly_degree}"](DOMAIN, [0., theta, delta_time])
  hdg.set_refinement(2 ** refinement)
  filename = f"leVeque_hyg{theta}-{poly_degree}-{refinement}"
  solution = common.time_stepping(
    hdg, time_steps, delta_time, solvers=("gmres", "bicgstab"),
    on_start=lambda x: common.plot(hdg, x, filename, time=time_end),
    on_step=lambda step, x: step % output_interval == 0 and hdg.plot_solution(list(x), time_end))
  error = hdg.errors(list(solution), time_end)[0]
  common.report(NAME, error, theta=theta, poly_degree=poly_degree, refinement=refinement)


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, theta=[0.5, 1.], poly_degree=range(3), refinement=range(5, 8))
