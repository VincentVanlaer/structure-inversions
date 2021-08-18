from argparse import ArgumentParser

code = "Gyre"
name = "simple-p-mode"
version = "1"

def apply_args(parser: ArgumentParser) -> None:
    parser.add_argument('--l', type=int, required=True)
    parser.add_argument('--m', type=int, required=True)

def get_suggested_model_name(args):
    return f'gyre_l{args.l}_m{args.m}_p'

def render(model, mesa_input, args) -> None:
    return f"""
&constants
/
&model
    model_type = 'EVOL'
    file = '{mesa_input}'
    file_format = 'MESA'
/

&rot
/

&mode
    l = {args.l}
    m = {args.m}
    tag = 'l1' ! Tag for namelist matching
    n_pg_min = 1
    n_pg_max = 60
/
&osc
    nonadiabatic = false
    outer_bound = 'UNNO'       ! Good choice as default
/
&num
    diff_scheme = 'COLLOC_GL4' ! second order scheme is good enough for non-adiabatic computations
    n_iter_max = 50
/

&scan
    grid_frame = 'INERTIAL'
    grid_type = 'LINEAR'      ! INVERSE for g-modes, LINEAR for p-modes
    freq_min = 0.1
    freq_max = 1000
    freq_max_units = 'CYC_PER_DAY'
    freq_min_units = 'CYC_PER_DAY'
    n_freq = 500               ! adjust as needed (increase if you notice missing radial orders in the list of consecutive modes)
    tag_list = 'l1'            ! Comma-separated list of tags to match
/

&grid
    w_osc = 10 ! Oscillatory region weight parameter  ! Good choice as default
    w_exp = 2  ! Exponential region weight parameter  ! Good choice as default
    w_ctr = 10 ! Central region weight parameter      ! Good choice as default
/


&ad_output
    freq_units = 'CYC_PER_DAY'
    freq_frame = 'INERTIAL'

    summary_file = '{model.joinpath('summary_ad_g_modes.HDF')}'     ! rename output file
    summary_file_format = 'HDF'
    summary_item_list = 'M_star,R_star,L_star,l,m,n_p,n_g,n_pg,omega,freq,E_norm,E,omega_int'    ! Items to appear in summary file (be careful not to put spaces in between items, this results in errors )

    detail_template = '{model.joinpath('ad_g_mode_l%L_n%N_j%J.HDF')}'
    detail_item_list = 'M_star,R_star,L_star,l,m,n,n_p,n_g,n_pg,freq,xi_h,xi_r,x,dW_dx,Gamma_1,P,rho,T,dE_dx,M_r,eul_P,eul_rho,deul_phi,eul_phi,dzeta_dx'

/

&nad_output
/
"""

def get_properties(args):
    return {
        'l': args.l,
        'm': args.m
    }
