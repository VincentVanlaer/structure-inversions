from argparse import ArgumentParser
#Gaël buldgen (Liege, genieva)
#Z, fov lower, envelope mixing
code = "Gyre"
name = "simple-g-mode"
version = "3"

def apply_args(parser: ArgumentParser) -> None:
    parser.add_argument('--l', type=int, required=True)
    parser.add_argument('--m', type=int, required=True)

def get_suggested_model_name(args):
    return f'l{args.l}_m{args.m}'

def render(model, mesa_input, args) -> str:
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
    n_pg_min = -60
    n_pg_max = -1
/
&osc
    nonadiabatic = false
    outer_bound = 'UNNO'       ! Good choice as default
/
&num
    diff_scheme = 'COLLOC_GL4'
    n_iter_max = 50
/

&scan
    grid_frame = 'INERTIAL'
    grid_type = 'INVERSE'      ! INVERSE for g-modes, LINEAR for p-modes
    freq_min = 0.1
    freq_max = 10
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
    detail_item_list = 'M_star,R_star,L_star,l,m,n,n_p,n_g,n_pg,freq,xi_h,xi_r,x,dW_dx,Gamma_1,P,rho,T,dE_dx,M_r,eul_P,eul_rho,deul_phi,eul_phi,dzeta_dx,As,c_1'

/

&nad_output
/
"""

def get_properties(args):
    return {
        'l': args.l,
        'm': args.m
    }
