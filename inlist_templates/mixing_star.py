from argparse import ArgumentParser
from pathlib import Path

code = "MESA"
name = "mixing-star"
version = "2"

def apply_args(parser: ArgumentParser) -> None:
    parser.add_argument('--mass', type=float, required=True)
    parser.add_argument('--overshoot', type=float, required=True)

def get_suggested_model_name(args) -> str:
    return f"{args.mass}M_fov_{args.overshoot}_Z002"

def render(model, args) -> str:
    if args.overshoot == 0:
        has_overshoot = '!'
    else:
        has_overshoot = ''

    return f"""
&star_job
  ! see star/defaults/star_job.defaults

  ! begin with a pre-main sequence model
    create_pre_main_sequence_model = .false.
    saved_model_name = '6M_002Z_at_ZAMS.mod'

  ! display on-screen plots
    pgstar_flag = .true.

    read_extra_star_job_inlist1 = .false.

    profile_columns_file = '{Path(__file__).parent.absolute().joinpath('profile_columns.list')}'

/ ! end of star_job namelist


&eos
  ! eos options
  ! see eos/defaults/eos.defaults
    read_extra_eos_inlist1 = .false.

/ ! end of eos namelist


&kap
  ! kap options
  ! see kap/defaults/kap.defaults
    use_Type2_opacities = .true.
    Zbase = 0.02

  ! Opacity tables to use
    kap_file_prefix = 'OP_a09_nans_removed_by_hand'       ! 'a09' vs 'OP_a09_nans_removed_by_hand'
    kap_CO_prefix = 'a09_co'
    kap_lowT_prefix = 'lowT_fa05_a09p'

/ ! end of kap namelist


&controls
  ! see star/defaults/controls.defaults
    scale_max_correction = 0.1 ! for pre-MS convergence

  ! starting specifications
    initial_mass = {args.mass} ! in Msun units
    initial_z = 0.02

  ! options for energy conservation (see MESA V, Section 3)
    use_dedt_form_of_energy_eqn = .true.
    use_gold_tolerances = .true.

    write_pulse_data_with_profile = .true.
    pulse_data_format = 'GYRE'
    log_directory = '{model.joinpath('LOGS')}'
    photo_directory = '{model.joinpath('photos')}'
    read_extra_controls_inlist1 = .false.

    use_Ledoux_criterion = .true.
    do_conv_premix = .true.
    num_cells_for_smooth_gradL_composition_term = 4
    num_cells_for_smooth_brunt_B = 4
    {has_overshoot}overshoot_D_min = 1d-2
    {has_overshoot}overshoot_f(1) = {args.overshoot}
    {has_overshoot}overshoot_f0(1) = 0.005
    {has_overshoot}overshoot_scheme(1) = 'exponential'
    {has_overshoot}overshoot_zone_type(1) = 'burn_H'
    {has_overshoot}overshoot_zone_loc = 'core'
    {has_overshoot}overshoot_bdy_loc = 'any'

  ! stop when the star nears ZAMS (Lnuc/L > 0.99)
    stop_near_zams = .false.

  ! stop when the center mass fraction of h1 drops below this limit
    xa_central_lower_limit_species(1) = 'h1'
    xa_central_lower_limit(1) = 1d-1

  ! Atmosphere
    atm_option = 'T_tau'
    atm_T_tau_relation = 'Eddington'
    atm_T_tau_opacity = 'varying'

  ! step & mesh control
    !time_delta_coeff = 0.15
    mesh_delta_coeff = 0.5
    max_allowed_nz = 20000
    okay_to_remesh = .true.
    use_other_mesh_delta_coeff_factor = .true. ! Overshoot regions
    use_other_mesh_functions = .true. ! Convective zones, set to false on the pre-MS since it can give issues there

    profile_interval = -1
    max_num_profile_models = -1
/ ! end of controls namelist

&pgstar

    read_extra_pgstar_inlist1 = .false.

/ ! end of pgstar namelist
"""

def get_properties(args):
    return {
        'mass': args.mass,
        'overshoot': args.overshoot
    }
