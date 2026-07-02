submodules/fiber_network.git/fiber_network/make_geo.py -i domains/fiber-2026-05-20/net1/sca -o output/test.geo.h5 --clamp .4
submodules/fiber_network.git/fiber_network/make_geo.py -i domains/fiber-2026-05-20/net1/sca -o output/test.geo --points output/test.pts --clamp .4

build/openblas/experiments/network -test constant -domain output/test.geo.h5  -ksp_norm_type unpreconditioned -ksp_monitor_yaml -mem_max  -pc_type net2as -net2as_pc_factor_mat_solver_type mumps -net2as_print_local -net2as_cb_trim -net2as_overlap_frac .1  -net2as_p 8 -net2as_cb_type q1
python experiments/ne18-08-dd-ref.py --domain output/test.geo --points output/test.pts --subdomains 8
