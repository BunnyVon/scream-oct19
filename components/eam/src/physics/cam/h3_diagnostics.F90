module h3_diagnostics

!-------------------------------------------------------------------------------
! Purpose:
! Gather h3 mass, surface fluxes, 
! then check global mean tritium mass conservation for different periods
! Conservation check frequency controlled by following namelist flags:
!   h3_print_diags_timestep - time step level checking (default on)
!   h3_print_diags_monthly  - monthly level checking   (default off)
!   h3_print_diags_total    - full run checking        (default on)
!
! Author: Sha Feng /20250423
!-------------------------------------------------------------------------------

use shr_kind_mod   , only: r8 => shr_kind_r8
use camsrfexch     , only: cam_in_t
use h3_cycle      , only:  h_i, h3_transport, h3_print_diags_timestep, &
                           h3_print_diags_monthly, h3_print_diags_total, &
                           h3_conserv_error_tol_per_year
use ppgrid         , only: pver, pcols, begchunk, endchunk
use physics_types  , only: physics_state, physics_tend, physics_ptend, &
                           physics_ptend_init
use constituents   , only: pcnst, cnst_name
use cam_logfile    , only: iulog
use spmd_utils     , only: masterproc
use cam_abortutils , only: endrun
use time_manager   , only: is_first_step, is_last_step, get_prev_date, &
                           get_curr_date, is_end_curr_month

implicit none
private
save

public h3_diags_init
public h3_diags_register
public get_total_h3
public get_h3_sfc_fluxes
public print_global_h3_diags
public h3_diags_store_fields
public h3_diags_read_fields

! Number of H3 tracers
integer, parameter :: ncnst = 2      ! number of H3 constituents

character(len=7), dimension(ncnst), parameter :: & ! constituent names
     c_names = (/'H3_ANT', 'H3    '/)

integer :: h3_ant_glo_ind ! global index of 'H3_ANT
integer :: h3_glo_ind     ! global index of 'H3'

!----- formats -----
character(*),parameter :: C_FA0   = "('    ',12x,(42x,a10,2x),' | ',(3x,a10,2x))"
character(*),parameter :: C_FF    = "('    ',a41,e25.17,' | ',e25.17)"
character(*),parameter :: C_FS_2  = "('    ',6x,e25.17,5x,e25.17,5x,' | ',e25.17)"
character(*),parameter :: C_FS2_2 = "('    ',a12,11x,e25.17,18x,' | ',e25.17)"
character(*),parameter :: C_SA0_2 = "('    ',8x,2(11x,a3,15x),' |',(8x,a12,2x))"
character(*),parameter :: C_RER   = "('    ',a15,8x,e25.17,21x)"


!-------------------------------------------------------------------------------
contains


   subroutine h3_diags_init(state)
      !-------------------------------------------------
      ! Purpose: initialize state h3 fields to zero
      !-------------------------------------------------

      !arguments
      type(physics_state), intent(inout) :: state(begchunk:endchunk)

      ! local variables
      integer :: ncol, lchnk, i

      ! zero out all of the state h3 variables
      do lchnk = begchunk, endchunk
         ncol = state(lchnk)%ncol
         do i = 1, ncol
            ! H3 at current, initial, month start, and previous time steps
            state(lchnk)%th_curr(i) = 0._r8
            state(lchnk)%th_init(i) = 0._r8
            state(lchnk)%th_mnst(i) = 0._r8
            state(lchnk)%th_prev(i) = 0._r8
            ! tritium emissions and fluxes at current timestep
            state(lchnk)%h_flux_sfc(i) = 0._r8
            state(lchnk)%h_flux_san(i) = 0._r8
            ! monthly accumulated tritium emissions and fluxes
            state(lchnk)%h_mflx_sfc(i) = 0._r8
            state(lchnk)%h_mflx_san(i) = 0._r8
            ! total time integrated tritium emissions and fluxes
            state(lchnk)%h_iflx_sfc(i) = 0._r8
            state(lchnk)%h_iflx_san(i) = 0._r8
         end do
      end do

   end subroutine h3_diags_init

!-------------------------------------------------------------------------------

   subroutine h3_diags_register()
      !-------------------------------------------------
      ! Purpose: register h3 fields into pbuf
      !-------------------------------------------------
      use physics_buffer, only: pbuf_add_field, dtype_r8

      integer :: idx

      if ( .not. h3_transport() ) return

      ! prior H3 amounts
      call pbuf_add_field('th_curr',   'global', dtype_r8, (/pcols/), idx)
      call pbuf_add_field('th_init',    'global', dtype_r8, (/pcols/), idx)
      call pbuf_add_field('th_mnst',    'global', dtype_r8, (/pcols/), idx)
      call pbuf_add_field('th_prev',    'global', dtype_r8, (/pcols/), idx)
      ! monthly accumulated h3 emissions and fluxes
      call pbuf_add_field('h_mflx_sfc', 'global', dtype_r8, (/pcols/), idx)
      call pbuf_add_field('h_mflx_san', 'global', dtype_r8, (/pcols/), idx)
      ! total accumulated h3 emissions and fluxes
      call pbuf_add_field('h_iflx_sfc', 'global', dtype_r8, (/pcols/), idx)
      call pbuf_add_field('h_iflx_san', 'global', dtype_r8, (/pcols/), idx)

   end subroutine h3_diags_register

!-------------------------------------------------------------------------------

   subroutine get_total_h3(state, wet_or_dry)
      !-------------------------------------------------
      ! Purpose: sum column h3 and store in state
      ! Called by: phys_run2
      !-------------------------------------------------
      use physconst,      only: rga

      type(physics_state), intent(inout) :: state
      character(len=3),    intent(in   ) :: wet_or_dry ! is h3 mmr wet or dry

      ! local variables
      real(r8) :: th(state%ncol)            ! vertical integral of total h3
      integer ncol                          ! number of atmospheric columns
      integer i, k, m                       ! column, level, constant indices
      !------------------------------------------------------------------------

      if ( .not. h3_transport() ) return

      ! Set H3 global index
      do m = 1, ncnst
         select case (trim(c_names(m)))
         case ('H3')
            h3_glo_ind = h_i(m)
         end select
      end do

      ! initialize array
      ncol = state%ncol
      do i = 1, ncol
         th(i) = 0._r8
      end do
      
      ! sum column h3 mass
      select case (trim(wet_or_dry))
      case ('wet')
         do k = 1, pver
            do i = 1, ncol
               th(i) = th(i) + state%q(i,k,h3_glo_ind) * state%pdel(i,k)
            end do
         end do
      case ('dry')
         do k = 1, pver
            do i = 1, ncol
               th(i) = th(i) + state%q(i,k,h3_glo_ind) * state%pdeldry(i,k)
            end do
         end do
      end select

      do i = 1, ncol
         th(i) = th(i) * rga
         print *, 'th(i) = ', i, th(i)
      end do

      ! copy new value to state
      do i = 1, ncol
         state%th_curr(i) = th(i)
      end do

   end subroutine get_total_h3

!-------------------------------------------------------------------------------

   subroutine get_h3_sfc_fluxes(state, cam_in, dtime)
      !-------------------------------------------------
      ! Purpose: store surface h3 exchange in state
      ! Called by: tphysac
      !-------------------------------------------------

      type(physics_state), intent(inout) :: state
      type(cam_in_t),      intent(in   ) :: cam_in
      real(r8), intent(in)               :: dtime        ! physics time step

      ! local variables
      real(r8) :: th(state%ncol)            ! vertical integral of total h3
      integer ncol                          ! number of atmospheric columns
      integer i, m                          ! column, constant indices
      real(r8) :: sfc_flux(pcols)           ! surface flux
      real(r8) :: sfc_flux_ant(pcols)       ! surface anthopogenic flux
      !------------------------------------------------------------------------

      if ( .not. h3_transport() ) return

      ! Set H3 global indices
      do m = 1, ncnst
         select case (trim(c_names(m)))
         case ('H3_ANT')
            h3_ant_glo_ind = h_i(m)
         case ('H3')
            h3_glo_ind     = h_i(m)
         end select
      end do

      ! initialize arrays
      ncol  = state%ncol
      do i = 1, ncol
         sfc_flux(i)     = 0._r8
         sfc_flux_ant(i) = 0._r8
      end do

      ! gather surface fluxes
      do i = 1, ncol
         sfc_flux(i)     = sfc_flux(i)     + cam_in%cflx(i,h3_glo_ind)
         sfc_flux_ant(i) = sfc_flux_ant(i) + cam_in%cflx(i,h3_ant_glo_ind)
      end do

      ! put in state
      do i = 1, ncol
         state%h_flux_sfc(i) = sfc_flux(i)
      end do

      ! zero out monthly fluxes at start of each month
      if ( is_start_curr_month() ) then
         do i = 1, ncol
            state%h_mflx_sfc(i) = 0._r8
            state%h_mflx_san(i) = 0._r8
         end do
      end if

      if ( .not. is_first_step() ) then
         do i = 1, ncol
            state%h_iflx_sfc(i) = state%h_iflx_sfc(i) + (sfc_flux(i)     * dtime)
            state%h_iflx_san(i) = state%h_iflx_san(i) + (sfc_flux_ant(i) * dtime)
            state%h_mflx_sfc(i) = state%h_mflx_sfc(i) + (sfc_flux(i)     * dtime)
            state%h_mflx_san(i) = state%h_mflx_san(i) + (sfc_flux_ant(i) * dtime)
         end do
      end if

   end subroutine get_h3_sfc_fluxes

!-------------------------------------------------------------------------------

   subroutine print_global_h3_diags(state, dtime, nstep)
      !-------------------------------------------------
      ! Purpose: Write out conservation checks
      ! Called by: phys_run2
      !-------------------------------------------------
      use phys_gmean,     only: gmean
      use phys_grid,      only: get_ncols_p

      type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
      real(r8), intent(in) :: dtime        ! physics time step
      integer , intent(in) :: nstep        ! current timestep number

      ! local variables
      integer :: ncol, lchnk, i
      integer :: ierr
      integer :: cdate, year, mon, day, sec
      integer, parameter :: h_num_var     = 4
      integer, parameter :: f_ts_num_var  = 1
      integer, parameter :: f_mon_num_var = 2
      integer, parameter :: f_run_num_var = 2
      character(len=*), parameter :: sub_name='print_global_h3_diags: '
      real(r8) :: time_integrated_flux, state_net_change
      real(r8) :: th_glob(h_num_var)
      real(r8) :: flux_ts_glob(f_ts_num_var)
      real(r8) :: flux_mon_glob(f_mon_num_var)
      real(r8) :: flux_run_glob(f_run_num_var)
      real(r8) :: gth_curr, gth_init, gth_mnst, gth_prev, gth_delta
      real(r8) :: gth_flux_sfc
      real(r8) :: gth_mflx_sfc, gth_mflx_san
      real(r8) :: gth_iflx_sfc, gth_iflx_san
      real(r8) :: gth_flux_tot, gth_mflx_tot, gth_iflx_tot
      real(r8) :: rel_error, expected_th
      real(r8) :: rel_error_mon, expected_th_mon
      real(r8) :: rel_error_run, expected_th_run
      real(r8) :: scaled_rel_error_tol, nyear
      real(r8), parameter :: seconds_per_year = 31536000._r8
      real(r8) :: seconds_in_month
      real(r8) :: total_seconds

      real(r8) :: th(      pcols,begchunk:endchunk,h_num_var)     ! array for holding h3 variables
      real(r8) :: flux_ts( pcols,begchunk:endchunk,f_ts_num_var)  ! array for holding timestep fluxes
      real(r8) :: flux_mon(pcols,begchunk:endchunk,f_mon_num_var) ! array for holding monthly fluxes
      real(r8) :: flux_run(pcols,begchunk:endchunk,f_run_num_var) ! array for holding full run fluxes
      !------------------------------------------------------------------------

      if ( .not. h3_transport() ) return

      do lchnk = begchunk, endchunk
         ncol = get_ncols_p(lchnk)
         do i = 1, ncol
            ! total h3 mass at different time points
            th(i,lchnk,1) = state(lchnk)%th_curr(i)
            th(i,lchnk,2) = state(lchnk)%th_init(i)
            th(i,lchnk,3) = state(lchnk)%th_mnst(i)
            th(i,lchnk,4) = state(lchnk)%th_prev(i)
            ! print *, 'th(i,lchnk) = ', i, lchnk, th(i,lchnk,1), th(i,lchnk,2), th(i,lchnk,3), th(i,lchnk,4)  
            ! h3 emissions and fluxes at current time step
            flux_ts(i,lchnk,1) = state(lchnk)%h_flux_sfc(i)
            ! flux_ts(i,lchnk,2) = state(lchnk)%h_flux_san(i)
            ! monthly accumulated h3 emissions and fluxes
            flux_mon(i,lchnk,1) = state(lchnk)%h_mflx_sfc(i)
            flux_mon(i,lchnk,2) = state(lchnk)%h_mflx_san(i)
            ! total time integrated h3 emissions and fluxes
            flux_run(i,lchnk,1) = state(lchnk)%h_iflx_sfc(i)
            flux_run(i,lchnk,2) = state(lchnk)%h_iflx_san(i)
         end do
      end do

      ! Compute global means of tritium variables
      if ( ( h3_print_diags_timestep ) .or. &
           ( h3_print_diags_monthly .and. is_end_curr_month() ) .or. & 
           ( h3_print_diags_total .and. is_last_step() ) ) then
         call gmean(th, th_glob, h_num_var)
         ! print *, 'th_glob = ', th_glob
      end if

      if ( h3_print_diags_timestep) then
         call gmean(flux_ts,  flux_ts_glob,  f_ts_num_var)
      end if
      ! print *, 'flux_ts_glob = ', flux_ts_glob
      if ( h3_print_diags_monthly .and. is_end_curr_month() ) then
         call gmean(flux_mon, flux_mon_glob, f_mon_num_var)
      end if
      ! print *, 'flux_mon_glob = ', flux_mon_glob
      if ( h3_print_diags_total .and. is_last_step() ) then
         call gmean(flux_run, flux_run_glob, f_run_num_var)
      end if
      ! print *, 'flux_run_glob = ', flux_run_glob

      ! assign global means to readable variables
      gth_curr     = th_glob(1)
      gth_init     = th_glob(2)
      gth_mnst     = th_glob(3)
      gth_prev     = th_glob(4)
      ! print *, 'gth_curr = ', gth_curr
      ! print *, 'gth_init = ', gth_init
      ! print *, 'gth_mnst = ', gth_mnst
      ! print *, 'gth_prev = ', gth_prev


      gth_flux_sfc = flux_ts_glob(1)

      gth_mflx_sfc = flux_mon_glob(1)
      gth_mflx_san = flux_mon_glob(2)

      gth_iflx_sfc = flux_run_glob(1)
      gth_iflx_san = flux_run_glob(2)

      ! Compute important terms
      gth_flux_tot    = gth_flux_sfc 
      gth_mflx_tot    = gth_mflx_sfc 
      gth_iflx_tot    = gth_iflx_sfc 
      ! print *, 'gth_flux_tot = ', gth_flux_tot
      
      expected_th     = gth_prev + (gth_flux_tot * dtime)
      ! write(iulog,*)'gth_curr = ',gth_curr
      ! write(iulog,*)'expected_th = ',expected_th
      
      rel_error       = ( expected_th - gth_curr ) / gth_curr

      expected_th_mon = gth_mnst + gth_mflx_tot ! dtime factor already included
      rel_error_mon   = (expected_th_mon - gth_curr) / gth_curr

      expected_th_run = gth_init + gth_iflx_tot ! dtime factor already included
      rel_error_run   = (expected_th_run - gth_curr) / gth_curr

      gth_delta    = gth_curr - gth_prev

      ! get the date
      call get_curr_date(year, mon, day, sec)
      cdate = year*10000 + mon*100 + day

      ! Time step level write outs----------------------------------------------
      if (masterproc .and. h3_print_diags_timestep) then

         write(iulog,*   )   ''
         write(iulog,*   )   'NET H3 FLUXES : period = timestep : date = ',cdate,sec
         write(iulog,C_FA0 ) '  Time  ',   '  Time    '
         write(iulog,C_FA0 ) 'averaged',   'integrated'
         write(iulog,C_FA0 ) 'kg/m2/s', 'kg/m2'

         write(iulog, '(71("-"),"|",20("-"))')

         write(iulog,C_FF) 'Surface  Emissions', gth_flux_sfc, gth_flux_sfc * dtime

         write(iulog, '(71("-"),"|",23("-"))')

         write(iulog,C_FF) '   *SUM*', &
              gth_flux_tot, gth_flux_tot * dtime

         time_integrated_flux = gth_flux_tot * dtime

         write(iulog, '(71("-"),"|",23("-"))')

         write(iulog,*)''
         write(iulog,*)'H3 MASS (kg/m2) : period = timestep : date = ',cdate,sec

         write(iulog,*)''
         write(iulog,C_SA0_2) 'beg', 'end', '*NET CHANGE*'
         write(iulog,C_FS_2) gth_prev, gth_curr, gth_delta


         write(iulog, '(71("-"),"|",23("-"))')

         write(iulog,C_FS2_2)'       *SUM*', &
              (gth_curr - gth_prev), &
              (gth_curr - gth_prev)

         state_net_change = (gth_curr - gth_prev)

         write(iulog,C_RER)'Relative Error:', rel_error

         if (nstep > 0) then
            if (abs(rel_error) > h3_conserv_error_tol_per_year) then
               write(iulog,*) 'time integrated flux = ', time_integrated_flux
               write(iulog,*) 'net change in state  = ', state_net_change
               write(iulog,*) 'error                = ', abs(time_integrated_flux - state_net_change)
               call endrun(trim(sub_name) // 'Mass conservation failure detected')
            end if
         end if

         write(iulog, '(71("-"),"|",23("-"))')
      end if ! (masterproc .and. h3_print_diags_timestep)

      ! Whole run write outs----------------------------------------------------
      if ( is_last_step() .and. h3_print_diags_total ) then
         total_seconds = nstep * dtime
         if (masterproc) then
            write(iulog,*   )   ''
            write(iulog,*   )   'NET H3 FLUXES : period = full run : date = ',cdate,sec
            write(iulog,C_FA0 ) '  Time  ',   '  Time    '
            write(iulog,C_FA0 ) 'averaged',   'integrated'
            write(iulog,C_FA0 ) 'kg/m2/s', 'kg/m2'

            write(iulog, '(71("-"),"|",20("-"))')

            write(iulog,C_FF) 'Accumulated Surface Flux      ', gth_iflx_sfc / total_seconds, gth_iflx_sfc

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FF) '   *SUM*', &
                 gth_iflx_tot / total_seconds, gth_iflx_tot

            time_integrated_flux = gth_iflx_tot

            write(iulog, '(71("-"),"|",20("-"))')
            write(iulog,C_FF) 'Accumulated Sfc H3 Ant Flux', gth_iflx_san / total_seconds, gth_iflx_san
            write(iulog, '(71("-"),"|",20("-"))')
            write(iulog,C_FF) '   *SUM*', &
                 (gth_iflx_san ) / total_seconds, &
                 (gth_iflx_san )

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,*)''
            write(iulog,*)'H3 MASS (kg/m2) : period = full run : date = ',cdate,sec

            write(iulog,*)''
            write(iulog,C_SA0_2) 'beg', 'end', '*NET CHANGE*'
            write(iulog,C_FS_2) gth_init, gth_curr, (gth_curr - gth_init)


            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FS2_2)'       *SUM*', &
                 (gth_curr - gth_init), &
                 (gth_curr - gth_init)

            state_net_change = (gth_curr - gth_init)

            write(iulog,C_RER)'Relative Error:', rel_error_run

            ! Allow error tolerance to grow in time for long simulation campaigns
            nyear = max(1._r8, (nstep * dtime) / seconds_per_year) ! set nyear to 1 during first year
            scaled_rel_error_tol = (1._r8 + h3_conserv_error_tol_per_year)**nyear - 1._r8

            if (nstep > 0) then
               if (abs(rel_error_run) > scaled_rel_error_tol) then
                  write(iulog,*) 'time integrated flux = ', time_integrated_flux
                  write(iulog,*) 'net change in state  = ', state_net_change
                  write(iulog,*) 'error                = ', abs(time_integrated_flux - state_net_change)
                  write(iulog,*) 'No point in erroring out now, but long-term h3 conservation is bad'
               end if
            end if

            write(iulog, '(71("-"),"|",23("-"))')
         end if ! (masterproc)
      end if ! ( is_last_step() .and. h3_print_diags_total )

      ! Monthly write outs------------------------------------------------------
      if ( is_end_curr_month() .and. h3_print_diags_monthly ) then
         call get_seconds_in_curr_month(seconds_in_month)
         if (masterproc) then
            write(iulog,*   )   ''
            write(iulog,*   )   'NET H3 FLUXES : period = monthly,: date = ',cdate,sec
            write(iulog,C_FA0 ) '  Time  ',   '  Time    '
            write(iulog,C_FA0 ) 'averaged',   'integrated'
            write(iulog,C_FA0 ) 'kg/m2/s', 'kg/m2'

            write(iulog, '(71("-"),"|",20("-"))')

            write(iulog,C_FF) 'Accumulated Surface Flux      ', gth_mflx_sfc / seconds_in_month, gth_mflx_sfc

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FF) '   *SUM*', &
                 gth_mflx_tot / seconds_in_month, gth_mflx_tot

            time_integrated_flux = gth_mflx_tot

            write(iulog, '(71("-"),"|",20("-"))')
            write(iulog,C_FF) 'Accumulated Sfc Ant Emiss', gth_mflx_san / seconds_in_month, gth_mflx_san
            write(iulog, '(71("-"),"|",20("-"))')
            write(iulog,C_FF) '   *SUM*', &
                 (gth_mflx_san ) / seconds_in_month, &
                 (gth_mflx_san )

            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,*)''
            write(iulog,*)'TRITIUM MASS (kg/m2) : period = monthly,: date = ',cdate,sec

            write(iulog,*)''
            write(iulog,C_SA0_2) 'beg', 'end', '*NET CHANGE*'
            write(iulog,C_FS_2) gth_mnst, gth_curr, (gth_curr - gth_mnst)


            write(iulog, '(71("-"),"|",23("-"))')

            write(iulog,C_FS2_2)'       *SUM*', &
                 (gth_curr - gth_mnst), &
                 (gth_curr - gth_mnst)

            state_net_change = (gth_curr - gth_mnst)

            write(iulog,C_RER)'Relative Error:', rel_error_mon

            if (nstep > 0) then
               if (abs(rel_error_mon) > h3_conserv_error_tol_per_year) then
                  write(iulog,*) 'time integrated flux = ', time_integrated_flux
                  write(iulog,*) 'net change in state  = ', state_net_change
                  write(iulog,*) 'error                = ', abs(time_integrated_flux - state_net_change)
                  call endrun(trim(sub_name) // 'Monthly conservation failure detected')
               end if
            end if

            write(iulog, '(71("-"),"|",23("-"))')
         end if ! (masterproc)
      end if ! ( is_end_curr_month() .and. h3_print_diags_monthly )

   end subroutine print_global_h3_diags

!-------------------------------------------------------------------------------

   subroutine h3_diags_store_fields(state, pbuf2d)
      !-------------------------------------------------
      ! Purpose: Store prior CO2 fields in physics buffer
      ! Called by: phys_run2
      !-------------------------------------------------
      use physics_types,  only: physics_state
      use ppgrid,         only: begchunk, endchunk
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index, pbuf_get_chunk


      !args
      type(physics_state), intent(in)        :: state(begchunk:endchunk)
      type(physics_buffer_desc), pointer     :: pbuf2d(:,:)

      !local vars
      type(physics_buffer_desc), pointer :: pbuf_chnk(:)

      integer  :: chnk, ncol, i
      real(r8), pointer, dimension(:) :: tmpptr_th_init
      real(r8), pointer, dimension(:) :: tmpptr_th_mnst
      real(r8), pointer, dimension(:) :: tmpptr_th_prev
      real(r8), pointer, dimension(:) :: tmpptr_h_mflx_sfc
      real(r8), pointer, dimension(:) :: tmpptr_h_mflx_san
      real(r8), pointer, dimension(:) :: tmpptr_h_iflx_sfc
      real(r8), pointer, dimension(:) :: tmpptr_h_iflx_san
      integer :: th_init_idx    = 0
      integer :: th_mnst_idx    = 0
      integer :: th_prev_idx    = 0
      integer :: h_mflx_sfc_idx = 0
      integer :: h_mflx_san_idx = 0
      integer :: h_iflx_sfc_idx = 0
      integer :: h_iflx_san_idx = 0

      if ( .not. h3_transport() ) return

      ! total H3
      th_init_idx    = pbuf_get_index('th_init')
      th_mnst_idx    = pbuf_get_index('th_mnst')
      th_prev_idx    = pbuf_get_index('th_prev')
      ! monthly fluxes
      h_mflx_sfc_idx = pbuf_get_index('h_mflx_sfc')
      h_mflx_san_idx = pbuf_get_index('h_mflx_san')
      ! total fluxes
      h_iflx_sfc_idx = pbuf_get_index('h_iflx_sfc')
      h_iflx_san_idx = pbuf_get_index('h_iflx_san')

      do chnk = begchunk,endchunk
         ncol = state(chnk)%ncol
         pbuf_chnk => pbuf_get_chunk(pbuf2d, chnk)
         ! total h3
         call pbuf_get_field(pbuf_chnk, th_init_idx, tmpptr_th_init )
         call pbuf_get_field(pbuf_chnk, th_mnst_idx, tmpptr_th_mnst )
         call pbuf_get_field(pbuf_chnk, th_prev_idx, tmpptr_th_prev )
         ! monthly fluxes
         call pbuf_get_field(pbuf_chnk, h_mflx_sfc_idx, tmpptr_h_mflx_sfc )
         call pbuf_get_field(pbuf_chnk, h_mflx_san_idx, tmpptr_h_mflx_san )
         ! total fluxes
         call pbuf_get_field(pbuf_chnk, h_iflx_sfc_idx, tmpptr_h_iflx_sfc )
         call pbuf_get_field(pbuf_chnk, h_iflx_san_idx, tmpptr_h_iflx_san )
         do i = 1, ncol
            ! total tritium
            tmpptr_th_init(i)    = state(chnk)%th_init(i)
            tmpptr_th_mnst(i)    = state(chnk)%th_mnst(i)
            tmpptr_th_prev(i)    = state(chnk)%th_prev(i)
            ! monthly fluxes
            tmpptr_h_mflx_sfc(i) = state(chnk)%h_mflx_sfc(i)
            tmpptr_h_mflx_san(i) = state(chnk)%h_mflx_san(i)
            ! total fluxes
            tmpptr_h_iflx_sfc(i) = state(chnk)%h_iflx_sfc(i)
            tmpptr_h_iflx_san(i) = state(chnk)%h_iflx_san(i)
         end do
      end do
   end subroutine h3_diags_store_fields

!-------------------------------------------------------------------------------

   subroutine h3_diags_read_fields(state, pbuf2d)
      !-------------------------------------------------
      ! Purpose: Retrieve prior CO2 fields and 
      !          set their appropriate state fields
      ! Called by: phys_run2
      !-------------------------------------------------
      use physics_types,  only: physics_state
      use ppgrid,         only: begchunk, endchunk
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index, pbuf_get_chunk

      type(physics_state), intent(inout) :: state(begchunk:endchunk)
      type(physics_buffer_desc), pointer :: pbuf2d(:,:)

      ! local variables
      type(physics_buffer_desc), pointer :: pbuf_chnk(:)
      integer ncol                           ! number of atmospheric columns
      integer chnk                           ! local chunk
      integer i                              ! column index
      real(r8), pointer, dimension(:) :: tmpptr_th_init
      real(r8), pointer, dimension(:) :: tmpptr_th_mnst
      real(r8), pointer, dimension(:) :: tmpptr_th_prev
      real(r8), pointer, dimension(:) :: tmpptr_h_mflx_sfc
      real(r8), pointer, dimension(:) :: tmpptr_h_mflx_san
      real(r8), pointer, dimension(:) :: tmpptr_h_iflx_sfc
      real(r8), pointer, dimension(:) :: tmpptr_h_iflx_san
      integer :: th_init_idx    = 0
      integer :: th_mnst_idx    = 0
      integer :: th_prev_idx    = 0
      integer :: h_mflx_sfc_idx = 0
      integer :: h_mflx_san_idx = 0
      integer :: h_iflx_sfc_idx = 0
      integer :: h_iflx_san_idx = 0
      !------------------------------------------------------------------------

      if ( .not. h3_transport() ) return

      ! acquire prior H3 totals from physics buffer
      th_init_idx    = pbuf_get_index('th_init')
      th_mnst_idx    = pbuf_get_index('th_mnst')
      th_prev_idx    = pbuf_get_index('th_prev')
      ! monthly fluxes
      h_mflx_sfc_idx = pbuf_get_index('h_mflx_sfc')
      h_mflx_san_idx = pbuf_get_index('h_mflx_san')
      ! total fluxes
      h_iflx_sfc_idx = pbuf_get_index('h_iflx_sfc')
      h_iflx_san_idx = pbuf_get_index('h_iflx_san')

      do chnk=begchunk,endchunk
         ncol = state(chnk)%ncol
         pbuf_chnk => pbuf_get_chunk(pbuf2d, chnk)
         ! total H3
         call pbuf_get_field(pbuf_chnk, th_init_idx, tmpptr_th_init )
         call pbuf_get_field(pbuf_chnk, th_mnst_idx, tmpptr_th_mnst )
         call pbuf_get_field(pbuf_chnk, th_prev_idx, tmpptr_th_prev )
         ! monthly fluxes
         call pbuf_get_field(pbuf_chnk, h_mflx_sfc_idx, tmpptr_h_mflx_sfc )
         call pbuf_get_field(pbuf_chnk, h_mflx_san_idx, tmpptr_h_mflx_san )
         ! total fluxes
         call pbuf_get_field(pbuf_chnk, h_iflx_sfc_idx, tmpptr_h_iflx_sfc )
         call pbuf_get_field(pbuf_chnk, h_iflx_san_idx, tmpptr_h_iflx_san )
         do i = 1, ncol
            ! total titrium
            state(chnk)%th_init(i)    = tmpptr_th_init(i)
            state(chnk)%th_mnst(i)    = tmpptr_th_mnst(i)
            state(chnk)%th_prev(i)    = tmpptr_th_prev(i)
            ! monthly fluxes
            state(chnk)%h_mflx_sfc(i) = tmpptr_h_mflx_sfc(i)
            state(chnk)%h_mflx_san(i) = tmpptr_h_mflx_san(i)
            ! total fluxes
            state(chnk)%h_iflx_sfc(i) = tmpptr_h_iflx_sfc(i)
            state(chnk)%h_iflx_san(i) = tmpptr_h_iflx_san(i)
         end do
      end do

   end subroutine h3_diags_read_fields

!-------------------------------------------------------------------------------

   logical function is_start_curr_month()
   ! Return true if current timestep is first of the current month
   ! Based on is_end_curr_month in time_manager.F90

   ! Local variables
     integer :: &
        yr,   &! year
        mon,  &! month
        day,  &! day of month
        tod    ! time of day (seconds past 00Z)

     call get_prev_date(yr, mon, day, tod)
     is_start_curr_month = (day == 1  .and.  tod == 0)
   end function is_start_curr_month


!-------------------------------------------------------------------------------

   subroutine get_seconds_in_curr_month(seconds_in_month)
   ! Return the number of seconds in the current month
   ! It is expected that this routine is
   ! called when is_end_curr_month is true

   ! Arguments
     real(r8), intent(out) :: seconds_in_month
   ! Local variables
     integer :: &
        yr,   &! year
        mon,  &! month
        day,  &! day of month
        tod    ! time of day (seconds past 00Z)
     real(r8), parameter :: seconds_per_day = 86400._r8

! if is_end_curr_month, then 
! get_prev_date should have day == last_day_of_month
     call get_prev_date(yr, mon, day, tod)
     seconds_in_month = seconds_per_day * day
     
   end subroutine get_seconds_in_curr_month

!-------------------------------------------------------------------------------

end module h3_diagnostics
