#*************************************************************
define constant epsilon = 0.075
define constant alpha = 20.0
define constant beta = 0.5
#*************************************************************
define constant pcfl_min = 1e20
define constant pcfl_max = 1e20
define constant pcfl_beta = 2
#*************************************************************
define constant newton_n = 10
define constant newton_res = 1.e-12

define constant min_order = 2
define constant max_order = 2

define constant stab_conv = 1
define constant stab_visc = 1
#*************************************************************
# PETSC
define constant show_monitor = 0
define constant ksp_restart = 120
define constant pc_factor_levels = 10
define constant rel_tol = 1.e-12
define constant max_steps = 1000
#*************************************************************
define constant output = none
define constant sol_write_error = 1
#*************************************************************
geometry = square.in2d
mesh = netgen_32.vol
#*************************************************************
define constant heapsize = 100000000
#*************************************************************
shared = "./main"
#*************************************************************
numproc Problem3and4 Problem3and4
