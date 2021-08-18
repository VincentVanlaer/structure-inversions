! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib

      implicit none

      ! Standard run_star_extras with addition of increased meshing by Siemen Burssens

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Begin insert by SB !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!

      ! declarations for extra meshing functions with inlist 'inlist_xtra_coeff_os_czb'
      real(dp) :: &
         xtra_coef_os_full_on, &
         xtra_coef_os_full_off, &
         xtra_coef_os_above_nonburn, &
         xtra_coef_os_below_nonburn, &
         xtra_coef_os_above_burn_h, &
         xtra_coef_os_below_burn_h, &
         xtra_coef_os_above_burn_he, &
         xtra_coef_os_below_burn_he, &
         xtra_coef_os_above_burn_z, &
         xtra_coef_os_below_burn_z, &
         xtra_dist_os_above_nonburn, &
         xtra_dist_os_below_nonburn, &
         xtra_dist_os_above_burn_h, &
         xtra_dist_os_below_burn_h, &
         xtra_dist_os_above_burn_he, &
         xtra_dist_os_below_burn_he, &
         xtra_dist_os_above_burn_z, &
         xtra_dist_os_below_burn_z, &
         xtra_mesh_czb_weight, &
         xtra_mesh_czb_width, &
         xtra_mesh_czb_center


      real(dp) :: Xc_save, Xc_save_step, Xc_precision
      integer :: Xc_count
      logical :: need_to_save

      namelist /xtra_coeff_os_czb/ &
         xtra_coef_os_full_on, &
         xtra_coef_os_full_off, &
         xtra_coef_os_above_nonburn, &
         xtra_coef_os_below_nonburn, &
         xtra_coef_os_above_burn_h, &
         xtra_coef_os_below_burn_h, &
         xtra_coef_os_above_burn_he, &
         xtra_coef_os_below_burn_he, &
         xtra_coef_os_above_burn_z, &
         xtra_coef_os_below_burn_z, &
         xtra_dist_os_above_nonburn, &
         xtra_dist_os_below_nonburn, &
         xtra_dist_os_above_burn_h, &
         xtra_dist_os_below_burn_h, &
         xtra_dist_os_above_burn_he, &
         xtra_dist_os_below_burn_he, &
         xtra_dist_os_above_burn_z, &
         xtra_dist_os_below_burn_z, &
         xtra_mesh_czb_weight, &
         xtra_mesh_czb_width, &
         xtra_mesh_czb_center

      ! end of declarations for xtra_coeff_os_czb

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! End insert by SB !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

    subroutine extras_controls(id, ierr)
       integer, intent(in) :: id
       integer, intent(out) :: ierr
       type (star_info), pointer :: s
       ierr = 0
       call star_ptr(id, s, ierr)
       if (ierr /= 0) return

       ! this is the place to set any procedure pointers you want to change
       ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


       ! the extras functions in this file will not be called
       ! unless you set their function pointers as done below.
       ! otherwise we use a null_ version which does nothing (except warn).

       !!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! Begin insert by SB !!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!
       Xc_save = 0.6       ! start saving models from this Xc value
       Xc_save_step = 0.005  ! interval in Xc to save
       Xc_precision = 1d-5
       Xc_count = 21

       s% extras_check_model => extras_check_model

       ! Meshing is now defined in run_star_extras through hook
       ! Note that two different hooks are called!
       ! OS is defined through other_mesh_delta_coeff_factor
       ! CZB is defined through other_mesh_fcn_dat

       call read_inlist_xtra_coeff_os_czb(ierr) ! Read inlist
       if (ierr /= 0) return
       s% how_many_other_mesh_fcns => how_many_mesh_fcns ! CZB
       s% other_mesh_fcn_data => gradr_grada_mesh_fcn_data ! CZB
       s% other_mesh_delta_coeff_factor => other_mesh_delta_coeff_factor ! CZB

       ! Once you have set the function pointers you want,
       ! then uncomment this (or set it in your star_job inlist)
       ! to disable the printed warning message,
       s% job% warn_run_star_extras = .false.

       !!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! End insert by SB !!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine extras_controls


    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Begin insert by SB !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine read_inlist_xtra_coeff_os_czb(ierr)
        use utils_lib
        integer, intent(out) :: ierr
        character (len=256) :: filename, message
        integer :: unit

        filename = 'inlist_xtra_coeff_os_czb'

        write(*,*) 'read_inlist_xtra_coeff_os_czb'

        ! set defaults
        xtra_coef_os_full_on = 1d-4
        xtra_coef_os_full_off = 0.1d0
        xtra_coef_os_above_nonburn = 1d0
        xtra_coef_os_below_nonburn = 1d0
        xtra_coef_os_above_burn_h = 1d0
        xtra_coef_os_below_burn_h = 1d0
        xtra_coef_os_above_burn_he = 1d0
        xtra_coef_os_below_burn_he = 1d0
        xtra_coef_os_above_burn_z = 1d0
        xtra_coef_os_below_burn_z = 1d0
        xtra_dist_os_above_nonburn = 0.2d0
        xtra_dist_os_below_nonburn = 0.2d0
        xtra_dist_os_above_burn_h = 0.2d0
        xtra_dist_os_below_burn_h = 0.2d0
        xtra_dist_os_above_burn_he = 0.2d0
        xtra_dist_os_below_burn_he = 0.2d0
        xtra_dist_os_above_burn_z = 0.2d0
        xtra_dist_os_below_burn_z = 0.2d0
        xtra_mesh_czb_weight = 500
        xtra_mesh_czb_width = 0.02d0
        xtra_mesh_czb_center = 0d0

        open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
        if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
        else
            read(unit, nml=xtra_coeff_os_czb, iostat=ierr)
            close(unit)
            if (ierr /= 0) then
                write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
                write(*, '(a)') &
                'The following runtime error message might help you find the problem'
                write(*, *)
                open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
                read(unit, nml=xtra_coeff_os_czb)
                close(unit)
            end if
        end if


    end subroutine read_inlist_xtra_coeff_os_czb

    ! Other mesh routine (copied from agb testsuite included with MESA v.12778)
    subroutine other_mesh_delta_coeff_factor(id, eps_h, eps_he, eps_z, ierr)
        use const_def
        use chem_def
        integer, intent(in) :: id
        real(dp), intent(in), dimension(:) :: eps_h, eps_he, eps_z
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        real(dp) :: he_cntr, full_off, full_on, alfa_os
        integer :: k, kk, nz, max_eps_loc
        real(dp) :: xtra_coef, xtra_dist, coef, Hp, r_extra, max_eps, eps
        logical :: in_convective_region
        logical, parameter :: dbg = .false.

        include 'formats'

        !write(*,*) 'enter other_mesh_delta_coeff_factor'
        ierr = 0
        if (xtra_coef_os_above_nonburn == 1d0 .and. &
        xtra_coef_os_below_nonburn == 1d0 .and. &
        xtra_coef_os_above_burn_h == 1d0 .and. &
        xtra_coef_os_below_burn_h == 1d0 .and. &
        xtra_coef_os_above_burn_he == 1d0 .and. &
        xtra_coef_os_below_burn_he == 1d0 .and. &
        xtra_coef_os_above_burn_z == 1d0 .and. &
        xtra_coef_os_below_burn_z == 1d0) return

        call star_ptr(id, s, ierr)
        if (ierr /= 0) return


        nz = s% nz
        he_cntr = s% xa(s% net_iso(ihe4),nz)
        full_off = xtra_coef_os_full_off
        full_on = xtra_coef_os_full_on
        if (he_cntr >= full_off) then
            alfa_os = 0
        else if (he_cntr <= full_on) then
            alfa_os = 1
        else
            alfa_os = (full_off - he_cntr)/(full_off - full_on)
        end if
        if (alfa_os == 0) return

        ! first go from surface to center doing below convective boundaries
        in_convective_region = (s% mixing_type(1) == convective_mixing)
        k = 2
        max_eps = -1d99
        max_eps_loc = -1
        do while (k <= nz)
            eps = eps_h(k) + eps_he(k) + eps_z(k)
            if (in_convective_region) then
                if (s% mixing_type(k) == convective_mixing) then
                    if (eps > max_eps) then
                        max_eps = eps
                        max_eps_loc = k
                    end if
                else
                    in_convective_region = .false.
                    if (max_eps < 1d0) then
                        xtra_coef = xtra_coef_os_below_nonburn
                        xtra_dist = xtra_dist_os_below_nonburn
                    else if (eps_h(max_eps_loc) > 0.5d0*max_eps) then
                        xtra_coef = xtra_coef_os_below_burn_h
                        xtra_dist = xtra_dist_os_below_burn_h
                    else if (eps_he(max_eps_loc) > 0.5d0*max_eps) then
                        xtra_coef = xtra_coef_os_below_burn_he
                        xtra_dist = xtra_dist_os_below_burn_he
                    else
                        xtra_coef = xtra_coef_os_below_burn_z
                        xtra_dist = xtra_dist_os_below_burn_z
                    end if
                    xtra_coef = xtra_coef*alfa_os + (1-alfa_os)
                    if (xtra_coef > 0 .and. xtra_coef /= 1) then
                        coef = xtra_coef
                        do
                            if (s% mixing_type(k) /= overshoot_mixing) exit
                            if (coef < s% mesh_delta_coeff_factor(k)) then
                                s% mesh_delta_coeff_factor(k) = coef
                            end if
                            if (k == nz) exit
                            k = k+1
                        end do
                        if (xtra_dist > 0) then
                            Hp = s% P(k)/(s% rho(k)*s% grav(k))
                            r_extra = max(0d0, s% r(k) - xtra_dist*Hp)
                            if (dbg) write(*,2) 'extra below overshoot region', &
                            k, s% r(k)/Rsun, Hp/Rsun, r_extra/Rsun
                            do
                                if (s% r(k) < r_extra) exit
                                if (coef < s% mesh_delta_coeff_factor(k)) then
                                    s% mesh_delta_coeff_factor(k) = coef
                                end if
                                if (k == nz) exit
                                k = k+1
                            end do
                        end if
                    end if
                    if (dbg) write(*,2) 'done with extra below overshoot region', k
                    if (dbg) write(*,*)
                end if
            else if (s% mixing_type(k) == convective_mixing) then
                in_convective_region = .true.
                max_eps = eps
                max_eps_loc = k
            end if
            k = k+1
        end do

        ! now go from center to surface doing above convective boundaries
        in_convective_region = (s% mixing_type(nz) == convective_mixing)
        k = nz-1
        max_eps = -1d99
        max_eps_loc = -1
        do while (k >= 1)
            eps = eps_h(k) + eps_he(k) + eps_z(k)
            if (in_convective_region) then
                if (s% mixing_type(k) == convective_mixing) then
                    if (eps > max_eps) then
                        max_eps = eps
                        max_eps_loc = k
                    end if
                else
                    in_convective_region = .false.
                    if (max_eps < 1d0) then
                        xtra_coef = xtra_coef_os_above_nonburn
                        xtra_dist = xtra_dist_os_above_nonburn
                    else if (eps_h(max_eps_loc) > 0.5d0*max_eps) then
                        xtra_coef = xtra_coef_os_above_burn_h
                        xtra_dist = xtra_dist_os_above_burn_h
                    else if (eps_he(max_eps_loc) > 0.5d0*max_eps) then
                        xtra_coef = xtra_coef_os_above_burn_he
                        xtra_dist = xtra_dist_os_above_burn_he
                    else
                        xtra_coef = xtra_coef_os_above_burn_z
                        xtra_dist = xtra_dist_os_above_burn_z
                    end if
                    xtra_coef = xtra_coef*alfa_os + (1-alfa_os)
                    if (dbg) write(*,2) 'xtra_coeff to surf', s% model_number, xtra_coef

                    if (xtra_coef > 0 .and. xtra_coef /= 1) then
                        coef = xtra_coef
                        do
                            if (s% mixing_type(k) /= overshoot_mixing) exit
                            if (coef < s% mesh_delta_coeff_factor(k)) then
                                s% mesh_delta_coeff_factor(k) = coef
                            end if
                            if (k == 1) exit
                            k = k-1
                        end do
                        if (xtra_dist > 0) then
                            Hp = s% P(k)/(s% rho(k)*s% grav(k))
                            r_extra = min(s% r(1), s% r(k) + xtra_dist*Hp)
                            if (dbg) write(*,2) 'extra above overshoot region', &
                            k, s% r(k)/Rsun, Hp/Rsun, r_extra/Rsun
                            do
                                if (s% r(k) > r_extra) exit
                                if (coef < s% mesh_delta_coeff_factor(k)) then
                                    s% mesh_delta_coeff_factor(k) = coef
                                end if
                                if (k == 1) exit
                                k = k-1
                            end do
                        end if
                    end if
                    if (dbg) write(*,2) 'done with extra above overshoot region', k
                    if (dbg) write(*,*)
                end if
            else if (s% mixing_type(k) == convective_mixing) then
                in_convective_region = .true.
                max_eps = eps
                max_eps_loc = k
            end if
            k = k-1
        end do

    end subroutine other_mesh_delta_coeff_factor

    subroutine how_many_mesh_fcns(id, n)
        integer, intent(in) :: id
        integer, intent(out) :: n
        n = 1
    end subroutine how_many_mesh_fcns

    subroutine gradr_grada_mesh_fcn_data( &
        id, nfcns, names, gval_is_xa_function, vals1, ierr)
        integer, intent(in) :: id
        integer, intent(in) :: nfcns
        character (len=*) :: names(:)
        logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
        real(dp), pointer :: vals1(:) ! =(nz, nfcns)
        integer, intent(out) :: ierr
        logical, parameter :: dbg = .false.

        real(dp), pointer :: vals(:,:)
        type (star_info), pointer :: s
        integer :: nz, k
        real(dp) :: weight, width, center
        real(dp), parameter :: maxval = 700d0 ! max value for tanh from crlibm

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        weight = xtra_mesh_czb_weight
        width =  xtra_mesh_czb_width
        center =  xtra_mesh_czb_center

        names(1) = 'gradr_grada_function'
        gval_is_xa_function(1) = .false.

        nz = s% nz
        vals(1:nz,1:nfcns) => vals1(1:nz*nfcns)

        ! Don't increase mesh when star settles on pre-MS (SB)
        ! Otherwise it blows up nz>40000
        if ((.not. s% doing_first_model_of_run) .and. (.not. s% doing_relax ) .and. (s% mixing_type(s%nz) .eq. convective_mixing)) then
            if (dbg) write(*,*) 'Mesh increased around czb'
            do k=1,nz
                vals(k,1) = weight*tanh(min(maxval,(s% gradr(k)-s% grada(k)-center)/width))*width
            end do
        end if

    end subroutine gradr_grada_mesh_fcn_data

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! End insert by SB !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    integer function extras_check_model(id)
        integer, intent(in) :: id
        integer :: ierr
        type (star_info), pointer :: s
        logical :: do_retry, do_terminate
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_check_model = keep_going

        do_retry = .false.
        do_terminate = .false.

        ! Terminate and save the pre-main sequence model when the convective core appears.
        ! The central hydrogen fraction needs to have decreased by a small amount,
        ! to make sure that core H-burning has started, and the star is near the ZAMS.
        ! If the model is not on the pre-main sequence, check if it needs to be saved,
        ! or if a retry needs to be made to save within precision at a desired Xc value.
        call save_at_Xc(id, Xc_save, Xc_precision, Xc_save_step, do_retry, do_terminate, ierr)
        if (do_retry) extras_check_model = retry
        if (do_terminate) extras_check_model = terminate

        ! by default, indicate where (in the code) MESA terminated
        if (extras_check_model == terminate) s% termination_code = t_extras_check_model
    end function extras_check_model

    subroutine save_at_Xc(id, Xc_save, Xc_precision, Xc_save_step, do_retry, do_terminate, ierr)
        integer, intent(in) :: id
        real(dp), intent(inout) :: Xc_save
        real(dp), intent(in) :: Xc_precision, Xc_save_step
        logical, intent(out) :: do_retry, do_terminate
        integer, intent(out) :: ierr
        type(star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (ierr /= 0) then
            write(*,*) 'Error: save_at_Xc: star_ptr failed'
            return
        endif

        do_retry  = .false.
        if (ABS(s% center_h1 - Xc_save) < Xc_precision) then
            Xc_count = Xc_count - 1
            if (Xc_count == 0) then
                do_terminate = .true.
            else
                Xc_save = Xc_save - Xc_save_step
                s% need_to_save_profiles_now = .true.
            endif
        else if (s% center_h1 < Xc_save) then
            s% dt = (s% xtra(1) - Xc_save) / (s% xtra(1) - s% center_h1) * s% dt
            do_retry = .true.
        endif
        s% xtra(1) = s% center_h1

    end subroutine save_at_Xc

end module run_star_extras
