#!/usr/bin/env python3
"""The parabolic advection solver (AdvectionParab, Parabolic loop, theta = 1/2, 10^3 time steps to
t = 10) with the injection test on the network domains/injection_test.geo, degree 3, with and
without the penalty tau: the state is plotted at the start and after every step into one file per
tau, the errors at t = 10 go to output/advection_hypergraph.yaml. Known since 2026-09: none
of the linear solvers converges, so shell_scripts/reproducibles_test.sh skips this script. The
classes AdvectionParab_D1_P3 come from modules/advection_injection_file_parabolic.cxx.

    python3 reproducibles_python/advection_hypergraph.py [--quick]
"""
import common

NAME = "advection_hypergraph"
hdg_module = common.module("advection_injection_file_parabolic")
DOMAIN = str(common.ROOT / "domains" / "injection_test.geo")


def diffusion_test(tau, theta, poly_degree, edge_dim, space_dim, iteration):
  time_steps, time_end = (10 ** 3 if not common.QUICK else 20), 10
  delta_time = time_end / time_steps
  hdg = hdg_module[f"AdvectionParab_D{edge_dim}_P{poly_degree}"](DOMAIN, [tau, theta, delta_time])
  filename = f"adv_injection{tau}"
  solution = common.time_stepping(
    hdg, time_steps, delta_time, solvers=("gmres", "bicgstab"),
    on_start=lambda x: common.plot(hdg, x, filename, time=time_end),
    on_step=lambda step, x: hdg.plot_solution(list(x), time_end))
  error = hdg.errors(list(solution), time_end)[0]
  common.report(NAME, error, tau=tau, theta=theta, poly_degree=poly_degree, edge_dim=edge_dim, space_dim=space_dim, iteration=iteration)


if __name__ == "__main__":
  common.main(diffusion_test, rtol=1e-5, tau=[0, 1], theta=[0.5], poly_degree=[3], edge_dim=[1], space_dim=[2], iteration=[2])
