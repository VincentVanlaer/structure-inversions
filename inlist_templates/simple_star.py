from argparse import ArgumentParser

code = "MESA"
name = "simple-star"
version = "1"

def apply_args(parser: ArgumentParser) -> None:
    parser.add_argument('--mass', type=float, required=True)

def render(model, args) -> None:
    return f"""
&star_job
  ! see star/defaults/star_job.defaults

  ! begin with a pre-main sequence model
    create_pre_main_sequence_model = .true.

  ! display on-screen plots
    pgstar_flag = .true.

    read_extra_star_job_inlist1 = .false.

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
    read_extra_kap_inlist1 = .false.

/ ! end of kap namelist


&controls
  ! see star/defaults/controls.defaults

  ! starting specifications
    initial_mass = {args.mass} ! in Msun units
    initial_z = 0.02

  ! options for energy conservation (see MESA V, Section 3)
     use_dedt_form_of_energy_eqn = .true.
     use_gold_tolerances = .true.

  ! stop when the center mass fraction of h1 drops below this limit
    xa_central_lower_limit_species(1) = 'h1'
    xa_central_lower_limit(1) = 1d-3

    write_pulse_data_with_profile = .true.
    pulse_data_format = 'GYRE'
    log_directory = '{model.joinpath('LOGS')}'
    photo_directory = '{model.joinpath('photos')}'
    read_extra_controls_inlist1 = .false.
/ ! end of controls namelist

&pgstar

    read_extra_pgstar_inlist1 = .false.

/ ! end of pgstar namelist
"""

def get_properties(args):
    return {
        'mass': args.mass
    }
