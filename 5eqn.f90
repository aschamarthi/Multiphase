! =============================================================================
! 5-equation compressible multiphase flow solver (2D)
!
! Equations solved: alpha*rho_1, alpha*rho_2, rho*u, rho*v, E, alpha_1
! Time integration: SSP-RK3 (Shu-Osher)
! Spatial flux:     HLLC Riemann solver with MP5 reconstruction (GRAB)
! Sainadh Chamarthi. Email: sainath@caltech.edu
! =============================================================================

program multiphase2d

    implicit none

    ! -------------------------------------------------------------------------
    ! Loop indices and counters
    ! -------------------------------------------------------------------------
    integer :: i, j, k
    integer :: N, rk_step
    integer :: NTMAX = 100000

    ! -------------------------------------------------------------------------
    ! Time variables
    ! -------------------------------------------------------------------------
    double precision :: t_end, time, dt
    double precision :: CFL = 0.4d0

    ! -------------------------------------------------------------------------
    ! Grid parameters
    ! -------------------------------------------------------------------------
    integer, parameter :: NX = 768*4, NY = 512*4, ghostp = 4, n_eqn = 6

    double precision :: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
    double precision :: dx, dy, xmin, xmax, ymin, ymax

    ! -------------------------------------------------------------------------
    ! Primitive variables
    ! -------------------------------------------------------------------------
    double precision :: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    ! -------------------------------------------------------------------------
    ! Phase variables: partial densities and volume fractions
    ! -------------------------------------------------------------------------
    double precision :: alpha_rho_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    ! -------------------------------------------------------------------------
    ! Equation of state parameters (stiffened gas)
    ! -------------------------------------------------------------------------
    double precision :: gamma_one, gamma_two, pi_one, pi_two
    double precision :: gamma_big, pi_big
    double precision :: gamma, p_inf

    ! -------------------------------------------------------------------------
    ! Conservative variables, old values, and residual
    ! -------------------------------------------------------------------------
    double precision :: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)
    double precision :: cons_old(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)
    double precision :: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)

    ! -------------------------------------------------------------------------
    ! Test case selector
    ! -------------------------------------------------------------------------
    integer :: test_case

    ! -------------------------------------------------------------------------
    ! Output controls
    ! -------------------------------------------------------------------------
    integer, parameter :: file_save = 600

    ! -------------------------------------------------------------------------
    ! Timing
    ! -------------------------------------------------------------------------
    integer          :: time_ini, time_end, time_calc
    double precision :: start, finish

    ! -------------------------------------------------------------------------
    ! Snapshot output times and flags
    ! -------------------------------------------------------------------------
    double precision :: t_end0, t_end1, t_end2, t_end3
    integer          :: flag_0, flag_1, flag_2, flag_3

    common /grid/ dx, dy

    ! =========================================================================
    ! Initialisation
    ! =========================================================================

    call cpu_time(start)
    write(*,*) 'Program start...'
    call system_clock(count = time_ini)

    flag_0 = 0;  flag_1 = 0;  flag_2 = 0;  flag_3 = 0

    ! Domain extents (in metres)
    xmin = 0.0d0
    xmax = 11.1d0/100.0d0

    ymin = 0.0d0
    ymax = 7.4d0/100.0d0

    ! -------------------------------------------------------------------------
    ! Build uniform Cartesian grid
    ! -------------------------------------------------------------------------
    dx = (xmax - xmin)/NX
    dy = (ymax - ymin)/NY

    do i = -ghostp, NX + ghostp
        x(i) = xmin + (i - 0.5d0)*dx
    enddo

    do j = -ghostp, NY + ghostp
        y(j) = ymin + (j - 0.5d0)*dy
    enddo

    ! =========================================================================
    ! Initial conditions
    ! =========================================================================

    call initialconditions(x, y, alpha_rho_one, alpha_rho_two, u_vel, v_vel,    &
                           pressure, alpha_one, NX, NY, test_case, t_end,        &
                           ghostp, n_eqn, gamma_one, gamma_two, pi_one, pi_two,  &
                           dx, dy)

    ! Convert primitive variables to conservatives
    do i = 1, NX
        do j = 1, NY

            alpha_two(i,j) = 1.0d0 - alpha_one(i,j)

            gamma_big = (alpha_one(i,j)/(gamma_one - 1.0d0))  &
                      + (alpha_two(i,j)/(gamma_two - 1.0d0))

            pi_big    = ((alpha_one(i,j)*gamma_one*pi_one)/(gamma_one - 1.0d0))  &
                      + ((alpha_two(i,j)*gamma_two*pi_two)/(gamma_two - 1.0d0))

            cons(i,j,1) = alpha_rho_one(i,j)
            cons(i,j,2) = alpha_rho_two(i,j)

            density(i,j) = alpha_rho_one(i,j) + alpha_rho_two(i,j)

            cons(i,j,3) = (alpha_rho_one(i,j) + alpha_rho_two(i,j))*u_vel(i,j)
            cons(i,j,4) = (alpha_rho_one(i,j) + alpha_rho_two(i,j))*v_vel(i,j)

            cons(i,j,5) = pressure(i,j)*gamma_big + pi_big  &
                        + 0.5*density(i,j)*(u_vel(i,j)**2 + v_vel(i,j)**2)

            cons(i,j,6) = alpha_one(i,j)

            gamma      = (1.0d0/gamma_big) + 1.0d0
            p_inf      = (gamma - 1.0d0)*pi_big/gamma
            sound(i,j) = (gamma*(pressure(i,j) + p_inf)/density(i,j))**0.5

        enddo
    enddo

    ! =========================================================================
    ! Time-loop setup
    ! =========================================================================

    time = 0.0d0

    ! Snapshot output times (seconds)
    t_end0 = 4d-6
    t_end1 = 17d-6
    t_end2 = 40d-6
    t_end3 = 67d-6

    N = 1
    call output(x, y, alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure,  &
                alpha_one, NX, NY, ghostp, N/file_save, time)

    call timestep(u_vel, v_vel, alpha_one, alpha_rho_one, alpha_rho_two,         &
                  density, gamma_one, gamma_two, pi_one, pi_two, pressure,        &
                  CFL, time, t_end, dt, NX, NY, ghostp, cons, n_eqn)

    call constoprim(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure,         &
                    alpha_one, sound, cons, NX, NY, ghostp, n_eqn,                &
                    gamma_one, gamma_two, pi_one, pi_two)

    write(*,*) dt, 'hey'

    write(*,*) '*********************************************'
    write(*,*) '   time step N        time             '
    write(*,*) '*********************************************'

    ! =========================================================================
    ! Main time integration loop  (SSP-RK3)
    ! =========================================================================

    do while (time .lt. t_end)

        cons_old = cons

        ! ----------------------------------------------------------------------
        ! RK stage 1:  cons^(1) = cons^n + dt * L(cons^n)
        ! ----------------------------------------------------------------------
        call boundaryconditions(alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                                pressure, alpha_one, x, y, time, NX, NY, ghostp)

        call FX(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure, alpha_one,  &
                residual, NX, NY, ghostp, n_eqn, dx, dy,                           &
                gamma_one, gamma_two, pi_one, pi_two)

        call GY(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure, alpha_one,  &
                residual, NX, NY, ghostp, n_eqn, dx, dy,                           &
                gamma_one, gamma_two, pi_one, pi_two)

        do k = 1, n_eqn
            do j = 1, NY
                do i = 1, NX
                    cons(i,j,k) = cons_old(i,j,k) + dt*residual(i,j,k)
                enddo
            enddo
        enddo

        call constoprim(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure,  &
                        alpha_one, sound, cons, NX, NY, ghostp, n_eqn,         &
                        gamma_one, gamma_two, pi_one, pi_two)

        ! ----------------------------------------------------------------------
        ! RK stage 2:  cons^(2) = 3/4 cons^n + 1/4 (cons^(1) + dt * L(cons^(1)))
        ! ----------------------------------------------------------------------
        call boundaryconditions(alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                                pressure, alpha_one, x, y, time, NX, NY, ghostp)

        call FX(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure, alpha_one,  &
                residual, NX, NY, ghostp, n_eqn, dx, dy,                           &
                gamma_one, gamma_two, pi_one, pi_two)

        call GY(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure, alpha_one,  &
                residual, NX, NY, ghostp, n_eqn, dx, dy,                           &
                gamma_one, gamma_two, pi_one, pi_two)

        do k = 1, n_eqn
            do j = 1, NY
                do i = 1, NX
                    cons(i,j,k) = (3.0d0/4.0d0)*cons_old(i,j,k)  &
                                + (1.0d0/4.0d0)*(cons(i,j,k) + dt*residual(i,j,k))
                enddo
            enddo
        enddo

        call constoprim(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure,  &
                        alpha_one, sound, cons, NX, NY, ghostp, n_eqn,         &
                        gamma_one, gamma_two, pi_one, pi_two)

        ! ----------------------------------------------------------------------
        ! RK stage 3:  cons^(n+1) = 1/3 cons^n + 2/3 (cons^(2) + dt * L(cons^(2)))
        ! ----------------------------------------------------------------------
        call boundaryconditions(alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                                pressure, alpha_one, x, y, time, NX, NY, ghostp)

        call FX(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure, alpha_one,  &
                residual, NX, NY, ghostp, n_eqn, dx, dy,                           &
                gamma_one, gamma_two, pi_one, pi_two)

        call GY(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure, alpha_one,  &
                residual, NX, NY, ghostp, n_eqn, dx, dy,                           &
                gamma_one, gamma_two, pi_one, pi_two)

        do k = 1, n_eqn
            do j = 1, NY
                do i = 1, NX
                    cons(i,j,k) = (1.0d0/3.0d0)*cons_old(i,j,k)  &
                                + (2.0d0/3.0d0)*(cons(i,j,k) + dt*residual(i,j,k))
                enddo
            enddo
        enddo

        call constoprim(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure,  &
                        alpha_one, sound, cons, NX, NY, ghostp, n_eqn,         &
                        gamma_one, gamma_two, pi_one, pi_two)

        ! ----------------------------------------------------------------------
        ! Advance time step
        ! ----------------------------------------------------------------------
        call timestep(u_vel, v_vel, alpha_one, alpha_rho_one, alpha_rho_two,    &
                      density, gamma_one, gamma_two, pi_one, pi_two, pressure,  &
                      CFL, time, t_end, dt, NX, NY, ghostp, cons, n_eqn)

        ! ----------------------------------------------------------------------
        ! Write snapshot outputs at requested times
        ! ----------------------------------------------------------------------
        if (flag_0 .eq. 1) then
            N = N + 800
            call output(x, y, alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                        pressure, alpha_one, NX, NY, ghostp, N/file_save, time)
            flag_0 = 2
            write(*,*) N/file_save, time
        endif

        if (flag_1 .eq. 1) then
            N = N + 800
            call output(x, y, alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                        pressure, alpha_one, NX, NY, ghostp, N/file_save, time)
            flag_1 = 2
            write(*,*) N/file_save, time
        endif

        if (flag_2 .eq. 1) then
            N = N + 800
            call output(x, y, alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                        pressure, alpha_one, NX, NY, ghostp, N/file_save, time)
            flag_2 = 2
            write(*,*) N/file_save, time
        endif

        if (flag_3 .eq. 1) then
            N = N + 800
            call output(x, y, alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                        pressure, alpha_one, NX, NY, ghostp, N/file_save, time)
            flag_3 = 2
            write(*,*) N/file_save, time
        endif

        ! Set flags and clip dt to hit snapshot times exactly
        if ((time+dt) .gt. t_end0 .and. flag_0 .eq. 0) then
            dt     = t_end0 - time
            flag_0 = 1
        endif

        if ((time+dt) .gt. t_end1 .and. flag_1 .eq. 0) then
            dt     = t_end1 - time
            flag_1 = 1
        endif

        if ((time+dt) .gt. t_end2 .and. flag_2 .eq. 0) then
            dt     = t_end2 - time
            flag_2 = 1
        endif

        if ((time+dt) .gt. t_end3 .and. flag_3 .eq. 0) then
            dt     = t_end3 - time
            flag_3 = 1
        endif

        if ((time+dt) .gt. t_end) dt = t_end - time

        time = time + dt
        N    = N + 1

        write(*,*) N, time, dt

    enddo

    ! =========================================================================
    ! Final output and timing
    ! =========================================================================

    write(*,*) '*********************************************'
    write(*,*) '   Number of time steps = ', N
    write(*,*) '*********************************************'

    call output(x, y, alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure,  &
                alpha_one, NX, NY, ghostp, N/file_save, time)

    call system_clock(count = time_end)
    time_calc = time_end - time_ini
    write(*,'(A20,I10,A)') 'Calculation time ', time_calc, ' [CPU ticks]'

    call cpu_time(finish)
    write(*,*) " Total CPU time to solution = ", finish-start, " seconds"
    write(*,*) 'Program ends...'

end program multiphase2d


! =============================================================================
! Inviscid flux in X-direction
!   Loops over all j-lines, reconstructs with GRAB (MP5), then applies
!   the HLLC Riemann solver at each i+1/2 face.
! =============================================================================

subroutine FX(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure, alpha_one,  &
              residual, NX, NY, ghostp, n_eqn, dx, dy,                           &
              gamma_one, gamma_two, pi_one, pi_two)

    implicit none

    integer          :: i, j, NX, NY, ghostp, k, n_eqn, ix
    double precision :: dx, dy

    double precision :: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision :: alpha_rho_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision :: gamma_one, gamma_two, pi_one, pi_two
    double precision :: gamma_big, pi_big
    double precision :: gamma, p_inf

    double precision :: primitive(-ghostp:NX+ghostp,n_eqn)
    double precision :: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)
    double precision :: source(-ghostp:NX+ghostp)

    double precision :: consl(-ghostp:NX+ghostp,n_eqn), consr(-ghostp:NX+ghostp,n_eqn)
    double precision :: prim_left(-ghostp:NX+ghostp,n_eqn), prim_right(-ghostp:NX+ghostp,n_eqn)

    double precision :: flux_half(-ghostp:NX+ghostp,n_eqn)
    double precision :: fright(-ghostp:NX+ghostp,n_eqn), fleft(-ghostp:NX+ghostp,n_eqn)

    double precision :: alpha_rho_one_left, alpha_rho_two_left, alpha_one_left
    double precision :: u_left, v_left, pressure_left, density_left, sound_left

    double precision :: alpha_rho_one_right, alpha_rho_two_right, alpha_one_right
    double precision :: u_right, v_right, pressure_right, density_right, sound_right

    double precision :: SL, SR, SP, EL, ER, u_avg, sound_avg
    double precision :: mx, my, minmod2

    double precision :: alpha_1, alpha_2, beta_1, beta_2
    double precision :: flux_hllc(-ghostp:NX+ghostp,n_eqn), flux_hll(-ghostp:NX+ghostp,n_eqn)

    double precision, parameter :: epsilon = 1.0d-5

    double precision :: indicator(-5:NX+5), cator2(-5:NX+5)
    double precision :: yuxin, aa, bb, sigma, f1, entropy(-5:NX+5)

    double precision :: blah = 0.35d0, crap = 1.0d-2

    ! Direction cosines: FX sweeps in x-direction
    mx = 1.0d0;  my = 0.0d0

    do j = 1, NY

        ! Pack 1D slice of primitive variables for reconstruction
        do i = -ghostp, NX + ghostp

            primitive(i,1) = alpha_rho_one(i,j)
            primitive(i,2) = alpha_rho_two(i,j)
            primitive(i,3) = u_vel(i,j)
            primitive(i,4) = v_vel(i,j)
            primitive(i,5) = pressure(i,j)
            primitive(i,6) = alpha_one(i,j)

            alpha_two(i,j) = 1.0d0 - alpha_one(i,j)

            gamma_big = ((alpha_one(i,j))/(gamma_one - 1.0d0))       &
                      + ((1.0d0 - alpha_one(i,j))/(gamma_two - 1.0d0))

            pi_big    = (((alpha_one(i,j))*gamma_one*pi_one)/(gamma_one - 1.0d0))          &
                      + (((1.0d0 - alpha_one(i,j))*gamma_two*pi_two)/(gamma_two - 1.0d0))

            gamma         = (1.0d0/gamma_big) + 1.0d0
            p_inf         = (gamma - 1.0d0)*pi_big/gamma
            density(i,j)  = alpha_rho_one(i,j) + alpha_rho_two(i,j)
            entropy(i)    = (pressure(i,j))/density(i,j)**gamma

        enddo

        ! Reconstruct left/right states at each i+1/2 face
        call GRAB(NX, primitive, prim_left, prim_right, mx, my, ghostp, n_eqn,  &
                  dx, gamma_one, gamma_two, pi_one, pi_two, cator2)

        ! ------------------------------------------------------------------
        ! HLLC Riemann solver at each x-face
        ! ------------------------------------------------------------------
        do i = -1, NX+1

            alpha_rho_one_left = prim_left(i,1)
            alpha_rho_two_left = prim_left(i,2)
            u_left             = prim_left(i,3)
            v_left             = prim_left(i,4)
            pressure_left      = prim_left(i,5)
            alpha_one_left     = prim_left(i,6)
            density_left       = prim_left(i,1) + prim_left(i,2)

            gamma_big = ((alpha_one_left)/(gamma_one - 1.0d0))              &
                      + ((1.0d0 - alpha_one_left)/(gamma_two - 1.0d0))

            pi_big    = (((alpha_one_left)*gamma_one*pi_one)/(gamma_one - 1.0d0))          &
                      + (((1.0d0 - alpha_one_left)*gamma_two*pi_two)/(gamma_two - 1.0d0))

            gamma = (1.0d0/gamma_big) + 1.0d0
            p_inf = (gamma - 1.0d0)*pi_big/gamma

            ! Conservative variables and fluxes on left state
            consl(i,1) = prim_left(i,1)
            consl(i,2) = prim_left(i,2)
            consl(i,3) = u_left*density_left
            consl(i,4) = v_left*density_left
            consl(i,5) = gamma_big*pressure_left + pi_big  &
                       + 0.5d0*(u_left**2.0d0 + v_left**2.0d0)*density_left
            consl(i,6) = prim_left(i,6)

            sound_left = (gamma*(pressure_left + p_inf)/density_left)**0.5

            fleft(i,1) = alpha_rho_one_left*u_left
            fleft(i,2) = alpha_rho_two_left*u_left
            fleft(i,3) = density_left*u_left**2.0d0 + pressure_left
            fleft(i,4) = density_left*u_left*v_left
            fleft(i,5) = (consl(i,5) + pressure_left)*u_left
            fleft(i,6) = alpha_one_left*u_left

            ! Right state
            alpha_rho_one_right = prim_right(i,1)
            alpha_rho_two_right = prim_right(i,2)
            u_right             = prim_right(i,3)
            v_right             = prim_right(i,4)
            pressure_right      = prim_right(i,5)
            alpha_one_right     = prim_right(i,6)
            density_right       = prim_right(i,1) + prim_right(i,2)

            gamma_big = ((alpha_one_right)/(gamma_one - 1.0d0))              &
                      + ((1.0d0 - alpha_one_right)/(gamma_two - 1.0d0))

            pi_big    = (((alpha_one_right)*gamma_one*pi_one)/(gamma_one - 1.0d0))          &
                      + (((1.0d0 - alpha_one_right)*gamma_two*pi_two)/(gamma_two - 1.0d0))

            gamma = (1.0d0/gamma_big) + 1.0d0
            p_inf = (gamma - 1.0d0)*pi_big/gamma

            consr(i,1) = prim_right(i,1)
            consr(i,2) = prim_right(i,2)
            consr(i,3) = u_right*density_right
            consr(i,4) = v_right*density_right
            consr(i,5) = gamma_big*pressure_right + pi_big  &
                       + 0.5d0*(u_right**2.0d0 + v_right**2.0d0)*density_right
            consr(i,6) = prim_right(i,6)

            sound_right = (gamma*(pressure_right + p_inf)/density_right)**0.5d0

            fright(i,1) = alpha_rho_one_right*u_right
            fright(i,2) = alpha_rho_two_right*u_right
            fright(i,3) = density_right*u_right**2.0d0 + pressure_right
            fright(i,4) = density_right*u_right*v_right
            fright(i,5) = (consr(i,5) + pressure_right)*u_right
            fright(i,6) = alpha_one_right*u_right

            ! Wave speed estimates (Einfeldt)
            u_avg     = (u_left + u_right)/2.0d0
            sound_avg = (sound_left + sound_right)/2.0d0

            SL = DMIN1(u_left  - sound_left,  u_avg - sound_avg)
            SR = DMAX1(u_right + sound_right, u_avg + sound_avg)

            SP = (pressure_right - pressure_left                                     &
               + density_left*u_left*(SL - u_left)                                   &
               - density_right*u_right*(SR - u_right))                               &
               / (density_left*(SL - u_left) - density_right*(SR - u_right))

            ! HLLC flux selection
            if (SL .gt. 0.0d0) then

                do k = 1, n_eqn
                    flux_hllc(i,k) = fleft(i,k)
                enddo

            else if (SL .le. 0.0d0 .and. 0.0d0 .lt. SP) then

                EL = (consl(i,5)) + (SP - u_left)*(SP*density_left + (pressure_left/(SL - u_left)))

                flux_hllc(i,1) = fleft(i,1) + SL*(prim_left(i,1)          *((SL-u_left)/(SL-SP)) - consl(i,1))
                flux_hllc(i,2) = fleft(i,2) + SL*(prim_left(i,2)          *((SL-u_left)/(SL-SP)) - consl(i,2))
                flux_hllc(i,3) = fleft(i,3) + SL*((density_left*SP)       *((SL-u_left)/(SL-SP)) - consl(i,3))
                flux_hllc(i,4) = fleft(i,4) + SL*((density_left*v_left)   *((SL-u_left)/(SL-SP)) - consl(i,4))
                flux_hllc(i,5) = fleft(i,5) + SL*(((SL-u_left)/(SL-SP))*EL                       - consl(i,5))
                flux_hllc(i,6) = fleft(i,6) + SL*(prim_left(i,6)          *((SL-u_left)/(SL-SP)) - consl(i,6))

            else if (SP .le. 0.0d0 .and. 0.0d0 .le. SR) then

                ER = (consr(i,5)) + (SP - u_right)*(SP*density_right + (pressure_right/(SR - u_right)))

                flux_hllc(i,1) = fright(i,1) + SR*(prim_right(i,1)        *((SR-u_right)/(SR-SP)) - consr(i,1))
                flux_hllc(i,2) = fright(i,2) + SR*(prim_right(i,2)        *((SR-u_right)/(SR-SP)) - consr(i,2))
                flux_hllc(i,3) = fright(i,3) + SR*((density_right*SP)     *((SR-u_right)/(SR-SP)) - consr(i,3))
                flux_hllc(i,4) = fright(i,4) + SR*((density_right*v_right)*((SR-u_right)/(SR-SP)) - consr(i,4))
                flux_hllc(i,5) = fright(i,5) + SR*(((SR-u_right)/(SR-SP))*ER                      - consr(i,5))
                flux_hllc(i,6) = fright(i,6) + SR*(prim_right(i,6)        *((SR-u_right)/(SR-SP)) - consr(i,6))

            elseif (SR .lt. 0.0d0) then

                do k = 1, n_eqn
                    flux_hllc(i,k) = fright(i,k)
                enddo

            end if

            flux_half(i,1:6) = flux_hllc(i,1:6)

            ! Non-conservative source term for volume fraction equation
            source(i) = ((1.0d0 + SIGN(1.0d0,SP))/2.0d0)*(u_left  + min(0.0d0,SL)*( ((SL-u_left )/(SL-SP)) - 1.0d0 )) &
                      + ((1.0d0 - SIGN(1.0d0,SP))/2.0d0)*(u_right + max(0.0d0,SR)*( ((SR-u_right)/(SR-SP)) - 1.0d0 ))

        enddo

        ! Divergence of flux -> residual
        do i = 1, NX
            residual(i,j,1:5) = -(flux_half(i,1:5) - flux_half(i-1,1:5))/dx
            residual(i,j,6)   = -(flux_half(i,6) - flux_half(i-1,6))/dx  &
                               + alpha_one(i,j)*(source(i) - source(i-1))/dx
        enddo

    enddo

end subroutine FX


! =============================================================================
! Inviscid flux in Y-direction
!   Mirrors FX but sweeps over i-lines and uses the v-velocity in the solver.
! =============================================================================

subroutine GY(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure, alpha_one,  &
              residual, NX, NY, ghostp, n_eqn, dx, dy,                           &
              gamma_one, gamma_two, pi_one, pi_two)

    implicit none

    integer          :: i, j, NX, NY, ghostp, k, n_eqn, ix
    double precision :: dx, dy

    double precision :: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision :: alpha_rho_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision :: gamma_one, gamma_two, pi_one, pi_two
    double precision :: gamma_big, pi_big
    double precision :: gamma, p_inf

    double precision :: primitive(-ghostp:NY+ghostp,n_eqn)
    double precision :: residual(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)
    double precision :: source(-ghostp:NY+ghostp)

    double precision :: consl(-ghostp:NY+ghostp,n_eqn), consr(-ghostp:NY+ghostp,n_eqn)
    double precision :: prim_left(-ghostp:NY+ghostp,n_eqn), prim_right(-ghostp:NY+ghostp,n_eqn)

    double precision :: flux_half(-ghostp:NY+ghostp,n_eqn)
    double precision :: fright(-ghostp:NY+ghostp,n_eqn), fleft(-ghostp:NY+ghostp,n_eqn)

    double precision :: alpha_rho_one_left, alpha_rho_two_left, alpha_one_left
    double precision :: u_left, v_left, pressure_left, density_left, sound_left

    double precision :: alpha_rho_one_right, alpha_rho_two_right, alpha_one_right
    double precision :: u_right, v_right, pressure_right, density_right, sound_right

    double precision :: SL, SR, SP, EL, ER, u_avg, sound_avg
    double precision :: mx, my, minmod2

    double precision :: alpha_1, alpha_2, beta_1, beta_2
    double precision :: flux_hllc(-ghostp:NY+ghostp,n_eqn), flux_hll(-ghostp:NY+ghostp,n_eqn)

    double precision, parameter :: epsilon = 1.0d-5

    double precision :: cator2(-5:NY+5), yuxin, aa, bb, sigma, f1, entropy(-5:NY+5)
    double precision :: blah = 0.35d0, crap = 1.0d-2

    ! Direction cosines: GY sweeps in y-direction
    mx = 0.0d0;  my = 1.0d0

    do i = 1, NX

        ! Pack 1D slice of primitive variables for reconstruction
        do j = -ghostp, NY + ghostp

            primitive(j,1) = alpha_rho_one(i,j)
            primitive(j,2) = alpha_rho_two(i,j)
            primitive(j,3) = u_vel(i,j)
            primitive(j,4) = v_vel(i,j)
            primitive(j,5) = pressure(i,j)
            primitive(j,6) = alpha_one(i,j)

            alpha_two(i,j) = 1.0d0 - alpha_one(i,j)

            gamma_big = ((alpha_one(i,j))/(gamma_one - 1.0d0))             &
                      + ((1.0d0 - alpha_one(i,j))/(gamma_two - 1.0d0))

            pi_big    = (((alpha_one(i,j))*gamma_one*pi_one)/(gamma_one - 1.0d0))          &
                      + (((1.0d0 - alpha_one(i,j))*gamma_two*pi_two)/(gamma_two - 1.0d0))

            gamma        = (1.0d0/gamma_big) + 1.0d0
            density(i,j) = alpha_rho_one(i,j) + alpha_rho_two(i,j)
            entropy(j)   = (pressure(i,j))/density(i,j)**gamma

        enddo

        ! Reconstruct left/right states at each j+1/2 face
        call GRAB(NY, primitive, prim_left, prim_right, mx, my, ghostp, n_eqn,  &
                  dy, gamma_one, gamma_two, pi_one, pi_two, cator2)

        ! ------------------------------------------------------------------
        ! HLLC Riemann solver at each y-face
        ! ------------------------------------------------------------------
        do j = -1, NY+1

            alpha_rho_one_left = prim_left(j,1)
            alpha_rho_two_left = prim_left(j,2)
            u_left             = prim_left(j,3)
            v_left             = prim_left(j,4)
            pressure_left      = prim_left(j,5)
            alpha_one_left     = prim_left(j,6)
            density_left       = prim_left(j,1) + prim_left(j,2)

            gamma_big = ((alpha_one_left)/(gamma_one - 1.0d0))              &
                      + ((1.0d0 - alpha_one_left)/(gamma_two - 1.0d0))

            pi_big    = (((alpha_one_left)*gamma_one*pi_one)/(gamma_one - 1.0d0))          &
                      + (((1.0d0 - alpha_one_left)*gamma_two*pi_two)/(gamma_two - 1.0d0))

            gamma = (1.0d0/gamma_big) + 1.0d0
            p_inf = (gamma - 1.0d0)*pi_big/gamma

            consl(j,1) = prim_left(j,1)
            consl(j,2) = prim_left(j,2)
            consl(j,3) = u_left*density_left
            consl(j,4) = v_left*density_left
            consl(j,5) = gamma_big*pressure_left + pi_big  &
                       + 0.5d0*(u_left**2.0d0 + v_left**2.0d0)*density_left
            consl(j,6) = prim_left(j,6)

            sound_left = (gamma*(pressure_left + p_inf)/density_left)**0.5

            ! Y-direction fluxes: momentum rows are transposed relative to FX
            fleft(j,1) = alpha_rho_one_left*v_left
            fleft(j,2) = alpha_rho_two_left*v_left
            fleft(j,3) = density_left*u_left*v_left
            fleft(j,4) = density_left*v_left**2.0d0 + pressure_left
            fleft(j,5) = (consl(j,5) + pressure_left)*v_left
            fleft(j,6) = alpha_one_left*v_left

            alpha_rho_one_right = prim_right(j,1)
            alpha_rho_two_right = prim_right(j,2)
            u_right             = prim_right(j,3)
            v_right             = prim_right(j,4)
            pressure_right      = prim_right(j,5)
            alpha_one_right     = prim_right(j,6)
            density_right       = prim_right(j,1) + prim_right(j,2)

            gamma_big = ((alpha_one_right)/(gamma_one - 1.0d0))              &
                      + ((1.0d0 - alpha_one_right)/(gamma_two - 1.0d0))

            pi_big    = (((alpha_one_right)*gamma_one*pi_one)/(gamma_one - 1.0d0))          &
                      + (((1.0d0 - alpha_one_right)*gamma_two*pi_two)/(gamma_two - 1.0d0))

            gamma = (1.0d0/gamma_big) + 1.0d0
            p_inf = (gamma - 1.0d0)*pi_big/gamma

            consr(j,1) = prim_right(j,1)
            consr(j,2) = prim_right(j,2)
            consr(j,3) = u_right*density_right
            consr(j,4) = v_right*density_right
            consr(j,5) = gamma_big*pressure_right + pi_big  &
                       + 0.5d0*(u_right**2.0d0 + v_right**2.0d0)*density_right
            consr(j,6) = prim_right(j,6)

            sound_right = (gamma*(pressure_right + p_inf)/density_right)**0.5

            fright(j,1) = alpha_rho_one_right*v_right
            fright(j,2) = alpha_rho_two_right*v_right
            fright(j,3) = density_right*u_right*v_right
            fright(j,4) = density_right*v_right**2.0d0 + pressure_right
            fright(j,5) = (consr(j,5) + pressure_right)*v_right
            fright(j,6) = alpha_one_right*v_right

            u_avg     = (v_left + v_right)/2.0d0
            sound_avg = (sound_left + sound_right)/2.0d0

            SL = DMIN1(v_left  - sound_left,  u_avg - sound_avg)
            SR = DMAX1(v_right + sound_right, u_avg + sound_avg)

            SP = (pressure_right - pressure_left                                      &
               + density_left*v_left*(SL - v_left)                                    &
               - density_right*v_right*(SR - v_right))                                &
               / (density_left*(SL - v_left) - density_right*(SR - v_right))

            ! HLLC flux selection
            if (SL .gt. 0.0d0) then

                do k = 1, n_eqn
                    flux_hllc(j,k) = fleft(j,k)
                enddo

            else if (SL .le. 0.0d0 .and. 0.0d0 .lt. SP) then

                EL = (consl(j,5)) + (SP - v_left)*(SP*density_left + (pressure_left/(SL - v_left)))

                flux_hllc(j,1) = fleft(j,1) + SL*(prim_left(j,1)         *((SL-v_left)/(SL-SP)) - consl(j,1))
                flux_hllc(j,2) = fleft(j,2) + SL*(prim_left(j,2)         *((SL-v_left)/(SL-SP)) - consl(j,2))
                flux_hllc(j,3) = fleft(j,3) + SL*((density_left*u_left)  *((SL-v_left)/(SL-SP)) - consl(j,3))
                flux_hllc(j,4) = fleft(j,4) + SL*((density_left*SP)      *((SL-v_left)/(SL-SP)) - consl(j,4))
                flux_hllc(j,5) = fleft(j,5) + SL*(((SL-v_left)/(SL-SP))*EL                      - consl(j,5))
                flux_hllc(j,6) = fleft(j,6) + SL*(prim_left(j,6)         *((SL-v_left)/(SL-SP)) - consl(j,6))

            else if (SP .le. 0.0d0 .and. 0.0d0 .le. SR) then

                ER = (consr(j,5)) + (SP - v_right)*(SP*density_right + (pressure_right/(SR - v_right)))

                flux_hllc(j,1) = fright(j,1) + SR*(prim_right(j,1)          *((SR-v_right)/(SR-SP)) - consr(j,1))
                flux_hllc(j,2) = fright(j,2) + SR*(prim_right(j,2)          *((SR-v_right)/(SR-SP)) - consr(j,2))
                flux_hllc(j,3) = fright(j,3) + SR*((density_right*u_right)  *((SR-v_right)/(SR-SP)) - consr(j,3))
                flux_hllc(j,4) = fright(j,4) + SR*((density_right*SP)       *((SR-v_right)/(SR-SP)) - consr(j,4))
                flux_hllc(j,5) = fright(j,5) + SR*(((SR-v_right)/(SR-SP))*ER                        - consr(j,5))
                flux_hllc(j,6) = fright(j,6) + SR*(prim_right(j,6)          *((SR-v_right)/(SR-SP)) - consr(j,6))

            elseif (SR .lt. 0.0d0) then

                do k = 1, n_eqn
                    flux_hllc(j,k) = fright(j,k)
                enddo

            end if

            flux_half(j,1:6) = flux_hllc(j,1:6)

            source(j) = ((1.0d0 + SIGN(1.0d0,SP))/2.0d0)*(v_left  + min(0.0d0,SL)*( ((SL-v_left )/(SL-SP)) - 1.0d0 )) &
                      + ((1.0d0 - SIGN(1.0d0,SP))/2.0d0)*(v_right + max(0.0d0,SR)*( ((SR-v_right)/(SR-SP)) - 1.0d0 ))

        enddo

        ! Accumulate y-direction contribution into residual
        do j = 1, NY
            residual(i,j,1:5) = residual(i,j,1:5) - (flux_half(j,1:5) - flux_half(j-1,1:5))/dy
            residual(i,j,6)   = residual(i,j,6)   - (flux_half(j,6) - flux_half(j-1,6))/dy  &
                               + alpha_one(i,j)*(source(j) - source(j-1))/dy
        enddo

    enddo

end subroutine GY


! =============================================================================
! GRAB - Wave-appropriate reconstruction subroutine
!   Reconstructs left/right interface states using MP5 in characteristic space.
!   Near liquid regions (p_inf > 2), falls back to MUSCL reconstruction.
!   Direction is controlled by (mx, my): (1,0) for x, (0,1) for y.
! =============================================================================

subroutine GRAB(NS, un, ulnew, urnew, mx, my, ghostp, n_eqn, ds,  &
                gamma_one, gamma_two, pi_one, pi_two, cator2)

    implicit none

    integer          :: ix, NS, ghostp, n_eqn, i, k, total
    double precision :: un(-ghostp:NS+ghostp,n_eqn), ds

    double precision :: ulnew(-ghostp:NS+ghostp,n_eqn)
    double precision :: urnew(-ghostp:NS+ghostp,n_eqn)
    double precision :: mx, my, lx, ly

    double precision :: lefteigen(n_eqn,n_eqn), righteigen(n_eqn,n_eqn)
    double precision :: v(-2:3,n_eqn), charstencil(-2:2), musclstencil(-1:2), vd2(-2:3,n_eqn)

    double precision, parameter :: epsilon = 1.0d-40, Constant = 1.0d00

    double precision :: ul(n_eqn), ur(n_eqn)

    ! MP5 constants
    double precision, PARAMETER :: B2 = 4.0d0/3.0d0, MP5 = 4.0d0, EPSM = 1.0d-20
    double precision :: VOR, VMP, DJM1, DJ, DJP1, DM4JPH, DM4JMH
    double precision :: VUL, VAV, VMD, VLC, VMIN, VMAX, MINMOD2, MINMOD4, l2norm, u_re

    double precision :: cator2(-5:NS+5)

    double precision :: vel, p, c, sqrt_rho_L, sqrt_rho_R, divisor, rho

    double precision :: alpha_rho_one_left, alpha_rho_one_right, alpha_rho_one_avg
    double precision :: alpha_rho_two_left, alpha_rho_two_right
    double precision :: alpha_rho_two_avg, alpha_one_avg, alpha_two_avg

    double precision :: gamma_one, gamma_two, pi_one, pi_two
    double precision :: gamma_big, pi_big
    double precision :: gamma, p_inf
    double precision :: del_m, del_o, del_p

    double precision, parameter :: kappa = 1.0d0/3.0d0

    double precision :: qa, qd, feng, T1, T2, CC, Beta_thinc
    double precision :: temp
    double precision :: wmin, wmax, wdelta, BB, AA, theta, avg

    ! Tangential direction cosines
    lx = -my;  ly = mx

    ! =========================================================================
    ! MP5 reconstruction in characteristic space
    ! =========================================================================

    do ix = 0, NS+1

        ! Compute Roe-averaged state at ix+1/2
        alpha_rho_one_left  = un(ix,1)
        alpha_rho_one_right = un(ix+1,1)
        alpha_rho_one_avg   = (alpha_rho_one_left + alpha_rho_one_right)/2.0d0

        alpha_rho_two_left  = un(ix,2)
        alpha_rho_two_right = un(ix+1,2)
        alpha_rho_two_avg   = (alpha_rho_two_left + alpha_rho_two_right)/2.0d0

        rho = alpha_rho_one_avg + alpha_rho_two_avg

        p = (un(ix,5) + un(ix+1,5))/2.0d0

        alpha_one_avg = (un(ix,6) + un(ix+1,6))/2.0d0
        alpha_two_avg = 1.0d0 - alpha_one_avg

        gamma_big = (alpha_one_avg/(gamma_one - 1.0d0))  &
                  + (alpha_two_avg/(gamma_two - 1.0d0))

        pi_big    = ((alpha_one_avg*gamma_one*pi_one)/(gamma_one - 1.0d0))  &
                  + ((alpha_two_avg*gamma_two*pi_two)/(gamma_two - 1.0d0))

        gamma = (1.0d0/gamma_big) + 1.0d0
        p_inf = (gamma - 1.0d0)*pi_big/gamma
        c     = (gamma*(p + p_inf)/rho)**0.5

        ! Left and right eigenvectors of the flux Jacobian
        righteigen(1,1) = alpha_rho_one_avg/(c*c*rho);  righteigen(1,2) = 1.0d0;  righteigen(1,3) = 0.0d0;  righteigen(1,4) = 0.0d0;  righteigen(1,5) = 0.0d0;  righteigen(1,6) = alpha_rho_one_avg/(c*c*rho)
        righteigen(2,1) = alpha_rho_two_avg/(c*c*rho);  righteigen(2,2) = 0.0d0;  righteigen(2,3) = 1.0d0;  righteigen(2,4) = 0.0d0;  righteigen(2,5) = 0.0d0;  righteigen(2,6) = alpha_rho_two_avg/(c*c*rho)
        righteigen(3,1) = -mx/(c*rho);                  righteigen(3,2) = 0.0d0;  righteigen(3,3) = 0.0d0;  righteigen(3,4) = my;     righteigen(3,5) = 0.0d0;  righteigen(3,6) = mx/(c*rho)
        righteigen(4,1) = -my/(c*rho);                  righteigen(4,2) = 0.0d0;  righteigen(4,3) = 0.0d0;  righteigen(4,4) = mx;     righteigen(4,5) = 0.0d0;  righteigen(4,6) = my/(c*rho)
        righteigen(5,1) = 1.0d0;                         righteigen(5,2) = 0.0d0;  righteigen(5,3) = 0.0d0;  righteigen(5,4) = 0.0d0;  righteigen(5,5) = 0.0d0;  righteigen(5,6) = 1.0d0
        righteigen(6,1) = 0.0d0;                         righteigen(6,2) = 0.0d0;  righteigen(6,3) = 0.0d0;  righteigen(6,4) = 0.0d0;  righteigen(6,5) = 1.0d0;  righteigen(6,6) = 0.0d0

        lefteigen(1,1) = 0.0d0;  lefteigen(1,2) = 0.0d0;  lefteigen(1,3) = -c*rho*mx*0.5d0;  lefteigen(1,4) = -c*rho*my*0.5d0;  lefteigen(1,5) = 0.5d0;                              lefteigen(1,6) = 0.0d0
        lefteigen(2,1) = 1.0d0;  lefteigen(2,2) = 0.0d0;  lefteigen(2,3) = 0.0d0;             lefteigen(2,4) = 0.0d0;            lefteigen(2,5) = -alpha_rho_one_avg/(c*c*rho);        lefteigen(2,6) = 0.0d0
        lefteigen(3,1) = 0.0d0;  lefteigen(3,2) = 1.0d0;  lefteigen(3,3) = 0.0d0;             lefteigen(3,4) = 0.0d0;            lefteigen(3,5) = -alpha_rho_two_avg/(c*c*rho);        lefteigen(3,6) = 0.0d0
        lefteigen(4,1) = 0.0d0;  lefteigen(4,2) = 0.0d0;  lefteigen(4,3) = my;                lefteigen(4,4) = mx;               lefteigen(4,5) = 0.0d0;                              lefteigen(4,6) = 0.0d0
        lefteigen(5,1) = 0.0d0;  lefteigen(5,2) = 0.0d0;  lefteigen(5,3) = 0.0d0;             lefteigen(5,4) = 0.0d0;            lefteigen(5,5) = 0.0d0;                              lefteigen(5,6) = 1.0d0
        lefteigen(6,1) = 0.0d0;  lefteigen(6,2) = 0.0d0;  lefteigen(6,3) = c*rho*mx*0.5d0;   lefteigen(6,4) = c*rho*my*0.5d0;  lefteigen(6,5) = 0.5d0;                              lefteigen(6,6) = 0.0d0

        ! ------------------------------------------------------------------
        ! Branch: MUSCL near liquid (p_inf > 2), MP5 otherwise
        ! ------------------------------------------------------------------
        if (p_inf .gt. 2.d0) then

            ! MUSCL with kappa = 1/3
            do i = -2, 3
                v(i,:) = un(i+ix,:)
            enddo

            do i = 1, n_eqn

                musclstencil = v(-1:2,i)

                del_m = musclstencil(0) - musclstencil(-1)
                del_o = musclstencil(1) - musclstencil(+0)
                del_p = musclstencil(2) - musclstencil(+1)

                ur(i) = musclstencil(1) - 0.25d0*( (1.0d0-kappa)*minmod2(del_p,2.0d0*del_o) &
                                                  + (1.0d0+kappa)*minmod2(del_o,2.0d0*del_p) )
                ul(i) = musclstencil(0) + 0.25d0*( (1.0d0-kappa)*minmod2(del_m,2.0d0*del_o) &
                                                  + (1.0d0+kappa)*minmod2(del_o,2.0d0*del_m) )

            enddo

            urnew(ix,:) = ur(:)
            ulnew(ix,:) = ul(:)

        else

            ! Project to characteristic variables
            do i = -2, 3
                v(i,:) = matmul(lefteigen,un(i+ix,:))
            enddo

            do i = 1, n_eqn

                ! Central 5-point average used for vortical wave (i==4)
                avg = ( ((1.0d0/60.0d0)*(2.0d0*v(-2,i) - 13.0d0*v(-1,i) + 47.0d0*v(0,i) + 27.0d0*v(1,i) - 3.0d0*v(2,i)))   &
                      + ((1.0d0/60.0d0)*(-3.0d0*v(-1,i) + 27.0d0*v(0,i) + 47.0d0*v(1,i) - 13.0d0*v(2,i) + 2.0d0*v(3,i))) ) &
                      * 0.5d0

                ! *** MP5 left reconstruction ***
                charstencil = v(-2:2,i)

                VOR = (1.0d0/60.0d0)*(2.0d0*charstencil(-2) - 13.0d0*charstencil(-1)  &
                    + 47.0d0*charstencil(0) + 27.0d0*charstencil(1) - 3.0d0*charstencil(2))

                if (i .eq. 4) VOR = avg

                VMP = charstencil(0) + MINMOD2(charstencil(1) - charstencil(0),  &
                                               MP5*(charstencil(0) - charstencil(-1)))

                if ((VOR - charstencil(0))*(VOR - VMP) < EPSM) then
                    ul(i) = VOR
                    if (i .eq. 4) ul(i) = avg
                else

                    DJM1 = charstencil(-2) - 2.0d0*charstencil(-1) + charstencil(0)
                    DJ   = charstencil(-1) - 2.0d0*charstencil(0)  + charstencil(1)
                    DJP1 = charstencil(0)  - 2.0d0*charstencil(1)  + charstencil(2)

                    DM4JPH = MINMOD4(4.0d0*DJ - DJP1, 4.0d0*DJP1 - DJ, DJ, DJP1)
                    DM4JMH = MINMOD4(4.0d0*DJ - DJM1, 4.0d0*DJM1 - DJ, DJ, DJM1)

                    VUL  = charstencil(0) + MP5*(charstencil(0) - charstencil(-1))
                    VAV  = 0.5d0*(charstencil(0) + charstencil(1))
                    VMD  = VAV  - 0.5d0*DM4JPH
                    VLC  = charstencil(0) + 0.5d0*(charstencil(0) - charstencil(-1)) + B2*DM4JMH

                    VMIN = MAX(MIN(charstencil(0),charstencil(1),VMD), MIN(charstencil(0),VUL,VLC))
                    VMAX = MIN(MAX(charstencil(0),charstencil(1),VMD), MAX(charstencil(0),VUL,VLC))

                    u_re = charstencil(0)  &
                         + (sign(1.0d0,charstencil(0) - charstencil(-1)) + sign(1.0d0,charstencil(1) - charstencil(0)))*0.5d0 &
                         * ((abs(charstencil(0)-charstencil(-1))*abs(charstencil(1)-charstencil(0))) &
                           /((abs(charstencil(0)-charstencil(-1)) + abs(charstencil(1)-charstencil(0))) + 1.0d-20))

                    ul(i) = (VOR + u_re)*0.5d0 - SIGN(1.0d0,(VOR-VMIN)*(VOR-VMAX))*(VOR - u_re)*0.5d0

                endif

                ! *** MP5 right reconstruction ***
                charstencil = v(-1:3,i)

                VOR = (1.0d0/60.0d0)*(-3.0d0*charstencil(-2) + 27.0d0*charstencil(-1)  &
                    + 47.0d0*charstencil(0) - 13.0d0*charstencil(1) + 2.0d0*charstencil(2))

                if (i .eq. 4) VOR = avg

                VMP = charstencil(0) + MINMOD2(charstencil(-1) - charstencil(0),  &
                                               MP5*(charstencil(0) - charstencil(1)))

                if ((VOR - charstencil(0))*(VOR - VMP) < EPSM) then
                    ur(i) = VOR
                    if (i .eq. 4) ur(i) = avg
                else

                    DJM1 = charstencil(-2) - 2.0d0*charstencil(-1) + charstencil(0)
                    DJ   = charstencil(-1) - 2.0d0*charstencil(0)  + charstencil(1)
                    DJP1 = charstencil(0)  - 2.0d0*charstencil(1)  + charstencil(2)

                    DM4JPH = MINMOD4(4.0d0*DJ - DJP1, 4.0d0*DJP1 - DJ, DJ, DJP1)
                    DM4JMH = MINMOD4(4.0d0*DJ - DJM1, 4.0d0*DJM1 - DJ, DJ, DJM1)

                    VUL  = charstencil(0) + MP5*(charstencil(0) - charstencil(+1))
                    VAV  = 0.5d0*(charstencil(0) + charstencil(-1))
                    VMD  = VAV  - 0.5d0*DM4JMH
                    VLC  = charstencil(0) + 0.5d0*(charstencil(0) - charstencil(+1)) + B2*DM4JPH

                    VMIN = MAX(MIN(charstencil(0),charstencil(-1),VMD), MIN(charstencil(0),VUL,VLC))
                    VMAX = MIN(MAX(charstencil(0),charstencil(-1),VMD), MAX(charstencil(0),VUL,VLC))

                    u_re = charstencil(0)  &
                         + (sign(1.0d0,-charstencil(0) + charstencil(-1)) + sign(1.0d0,-charstencil(1) + charstencil(0)))*0.5d0 &
                         * ((abs(charstencil(0)-charstencil(-1))*abs(charstencil(1)-charstencil(0))) &
                           /((abs(charstencil(0)-charstencil(-1)) + abs(charstencil(1)-charstencil(0))) + 1.0d-20))

                    ur(i) = (VOR + u_re)*0.5d0 - SIGN(1.0d0,(VOR-VMIN)*(VOR-VMAX))*(VOR - u_re)*0.5d0

                endif

            enddo

            ! Project back to physical space
            urnew(ix,:) = matmul(righteigen,ur)
            ulnew(ix,:) = matmul(righteigen,ul)

        endif

        ! Enforce positivity of partial densities and volume fraction
        if (urnew(ix,1) .lt. 0.0d0) urnew(ix,1) = un(ix+1,1)
        if (urnew(ix,2) .lt. 0.0d0) urnew(ix,2) = un(ix+1,2)
        if (urnew(ix,6) .lt. 0.0d0 .or. urnew(ix,6) .gt. 1.0d0) urnew(ix,6) = un(ix+1,6)

        if (ulnew(ix,1) .lt. 0.0d0) ulnew(ix,1) = un(ix,1)
        if (ulnew(ix,2) .lt. 0.0d0) ulnew(ix,2) = un(ix,2)
        if (ulnew(ix,6) .lt. 0.0d0 .or. ulnew(ix,6) .gt. 1.0d0) ulnew(ix,6) = un(ix,6)

        ! Revert to cell averages if speed of sound is imaginary
        if (c .lt. 0.0d0) then
            urnew(ix,:) = un(ix+1,:)
            ulnew(ix,:) = un(ix,:)
        endif

    enddo

end subroutine GRAB


! =============================================================================
! minmod2 - two-argument minmod limiter (used by MP5)
! =============================================================================

double precision function minmod2(X,Y)
    double precision, INTENT(IN) :: X, Y
    minmod2 = 0.5d0*(sign(1.0d0,X) + sign(1.0d0,Y))*min(ABS(X),ABS(Y))
end function minmod2


! =============================================================================
! minmod4 - four-argument minmod limiter (used by MP5)
! =============================================================================

double precision function minmod4(W,X,Y,Z)
    double precision, INTENT(IN) :: W, X, Y, Z

    minmod4 = 0.125d0*(sign(1.0d0,W) + sign(1.0d0,X))                       &
            * ABS((sign(1.0d0,W) + sign(1.0d0,Y))*(sign(1.0d0,W) + SIGN(1.0d0,Z))) &
            * min(ABS(W),ABS(X),ABS(Y),ABS(Z))

end function minmod4


! =============================================================================
! Boundary conditions
!   Left/right: outflow (zero-gradient)
!   Top/bottom: reflecting wall (v_vel negated)
! =============================================================================

subroutine boundaryconditions(alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                               pressure, alpha_one, x, y, time, NX, NY, ghostp)

    implicit none

    integer          :: i, j, NX, NY, ghostp
    double precision :: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp)
    double precision :: dx, dy
    double precision :: time
    double precision :: pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision, parameter :: pi = acos(-1.0d0)

    common /grid/ dx, dy

    ! Left and right ghost cells (zero-gradient / outflow)
    do i = 1, ghostp
        do j = -ghostp, NY+ghostp

            alpha_rho_one(-i+1,j)  = alpha_rho_one(1,j)
            alpha_rho_two(-i+1,j)  = alpha_rho_two(1,j)
            u_vel(-i+1,j)          = u_vel(1,j)
            v_vel(-i+1,j)          = v_vel(1,j)
            pressure(-i+1,j)       = pressure(1,j)
            alpha_one(-i+1,j)      = alpha_one(1,j)

            alpha_rho_one(NX+i,j)  = alpha_rho_one(NX,j)
            alpha_rho_two(NX+i,j)  = alpha_rho_two(NX,j)
            u_vel(NX+i,j)          = u_vel(NX,j)
            v_vel(NX+i,j)          = v_vel(NX,j)
            pressure(NX+i,j)       = pressure(NX,j)
            alpha_one(NX+i,j)      = alpha_one(NX,j)

        enddo
    enddo

    ! Bottom and top ghost cells (reflecting wall)
    do i = -ghostp, NX+ghostp
        do j = 1, ghostp

            alpha_rho_one(i,-j+1)  = alpha_rho_one(i,1)
            alpha_rho_two(i,-j+1)  = alpha_rho_two(i,1)
            u_vel(i,-j+1)          = u_vel(i,1)
            v_vel(i,-j+1)          = -v_vel(i,1)
            pressure(i,-j+1)       = pressure(i,1)
            alpha_one(i,-j+1)      = alpha_one(i,1)

            alpha_rho_one(i,NY+j)  = alpha_rho_one(i,NY)
            alpha_rho_two(i,NY+j)  = alpha_rho_two(i,NY)
            u_vel(i,NY+j)          = u_vel(i,NY)
            v_vel(i,NY+j)          = -v_vel(i,NY)
            pressure(i,NY+j)       = pressure(i,NY)
            alpha_one(i,NY+j)      = alpha_one(i,NY)

        enddo
    enddo

end subroutine boundaryconditions


! =============================================================================
! Conservative to primitive variable conversion
!   Recovers alpha*rho_k, u, v, p, alpha_1 and speed of sound from cons.
! =============================================================================

subroutine constoprim(alpha_rho_one, alpha_rho_two, u_vel, v_vel, pressure,  &
                      alpha_one, sound, cons, NX, NY, ghostp, n_eqn,         &
                      gamma_one, gamma_two, pi_one, pi_two)

    implicit none

    integer :: i, j, ghostp, NX, NY, n_eqn

    double precision :: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision :: alpha_rho_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision :: gamma_one, gamma_two, pi_one, pi_two
    double precision :: gamma_big, pi_big
    double precision :: gamma, p_inf

    double precision :: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)

    do i = 1, NX
        do j = 1, NY

            alpha_rho_one(i,j) = cons(i,j,1)
            alpha_rho_two(i,j) = cons(i,j,2)

            alpha_one(i,j) = cons(i,j,6)

            ! Clip volume fraction to [epsilon, 1-epsilon]
            if (alpha_one(i,j) .lt. 0.0) alpha_one(i,j) = 1.0d-8
            if (alpha_one(i,j) .gt. 1.0) alpha_one(i,j) = (1.0d00 - 1.0d-8)

            alpha_two(i,j) = 1.0d0 - alpha_one(i,j)

            gamma_big = (alpha_one(i,j)/(gamma_one - 1.0d0))  &
                      + (alpha_two(i,j)/(gamma_two - 1.0d0))

            pi_big    = ((alpha_one(i,j)*gamma_one*pi_one)/(gamma_one - 1.0d0))  &
                      + ((alpha_two(i,j)*gamma_two*pi_two)/(gamma_two - 1.0d0))

            density(i,j) = cons(i,j,1) + cons(i,j,2)

            u_vel(i,j)    = cons(i,j,3)/density(i,j)
            v_vel(i,j)    = cons(i,j,4)/density(i,j)

            pressure(i,j) = (cons(i,j,5) - pi_big  &
                           - 0.5*density(i,j)*(u_vel(i,j)**2 + v_vel(i,j)**2)) / gamma_big

            gamma      = (1.0d0/gamma_big) + 1.0d0
            p_inf      = (gamma - 1.0d0)*pi_big/gamma
            sound(i,j) = (gamma*(pressure(i,j) + p_inf)/density(i,j))**0.5

        enddo
    enddo

end subroutine constoprim


! =============================================================================
! Tridiagonal solver (Thomas algorithm)
!   Solves l*x(i-1) + d*x(i) + cons*x(i+1) = r(i) for x
! =============================================================================

subroutine tridiag(l, d, cons, r, x, n)

    implicit none

    integer,                        intent(in)  :: n
    double precision, dimension(n), intent(in)  :: l, d, cons, r
    double precision, dimension(n), intent(out) :: x

    integer          :: i, j
    double precision :: t
    double precision, dimension(n) :: g

    ! Forward sweep (LU decomposition + substitution)
    t    = d(1)
    x(1) = r(1)/t

    do i = 2, n
        j    = i - 1
        g(i) = cons(j)/t
        t    = d(i) - l(i)*g(i)

        if (t == 0.0d0) then
            write(*,*) "algebra::tridiag", "solution failed!"
            stop
            ! Note: with ifort this singularity does not occur; gfortran bug workaround:
            ! t = 1.0d-40
        end if

        x(i) = (r(i) - l(i)*x(j))/t
    end do

    ! Backward substitution
    do i = n - 1, 1, -1
        j    = i + 1
        x(i) = x(i) - g(j)*x(j)
    end do

end subroutine tridiag


! =============================================================================
! Time-step calculation (CFL condition)
!   Computes the maximum stable dt based on acoustic CFL.
! =============================================================================

subroutine timestep(u_vel, v_vel, alpha_one, alpha_rho_one, alpha_rho_two,      &
                    density, gamma_one, gamma_two, pi_one, pi_two, pressure,     &
                    CFL, time, t_end, dt, NX, NY, ghostp, cons, n_eqn)

    implicit none

    integer          :: i, j
    integer          :: NX, NY, ghostp, n_eqn
    double precision :: dx, dy
    double precision :: dt, time, t_end, dtnew, CFL

    double precision :: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: cons(-ghostp:NX+ghostp,-ghostp:NY+ghostp,n_eqn)

    double precision :: gamma_one, gamma_two, pi_one, pi_two
    double precision :: gamma_big, pi_big
    double precision :: gamma, p_inf

    double precision :: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: x_velocity, y_velocity

    double precision :: mu_lam, dt_visc

    common /grid/ dx, dy

    mu_lam  = 1.0d0/3.0d5
    dt_visc = 0.25d0*DMIN1(dx**2/mu_lam, dy**2/mu_lam)
    dt      = 1.0d10

    x_velocity = 0.0d0
    y_velocity = 0.0d0

    do i = 1, NX
        do j = 1, NY

            alpha_two(i,j) = 1.0d0 - alpha_one(i,j)

            gamma_big = (alpha_one(i,j)/(gamma_one - 1.0d0))  &
                      + (alpha_two(i,j)/(gamma_two - 1.0d0))

            pi_big    = ((alpha_one(i,j)*gamma_one*pi_one)/(gamma_one - 1.0d0))  &
                      + ((alpha_two(i,j)*gamma_two*pi_two)/(gamma_two - 1.0d0))

            density(i,j) = alpha_rho_one(i,j) + alpha_rho_two(i,j)

            gamma      = (1.0d0/gamma_big) + 1.0d0
            p_inf      = (gamma - 1.0d0)*pi_big/gamma
            sound(i,j) = (gamma*(pressure(i,j) + p_inf)/density(i,j))**0.5

            x_velocity = dx/max(x_velocity, ABS(u_vel(i,j)) + sound(i,j))
            y_velocity = dy/max(y_velocity, ABS(v_vel(i,j)) + sound(i,j))

            ! Multidimensional CFL combination
            dtnew = x_velocity*y_velocity/(x_velocity + y_velocity)
            if (dtnew .lt. dt) dt = dtnew

        enddo
    enddo

    dt = CFL*dt

    ! Do not overshoot end time
    if ((time + dt) .gt. t_end) dt = t_end - time

end subroutine timestep


! =============================================================================
! Output subroutine
!   Writes Tecplot .plt files with density, velocity, pressure, volume
!   fraction and a Schlieren-like density gradient magnitude field.
! =============================================================================

subroutine output(x, y, alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                  pressure, alpha_one, NX, NY, ghostp, file, time)

    implicit none

    integer          :: i, j, file
    integer          :: NX, NY, ghostp
    double precision :: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp), time, dx, dy

    double precision :: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    character(len=8) :: number*4, file_name

    double precision :: density_gradx(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: phi(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: density_grady(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: dg_max, density_grad_mag(1:NX,1:NY)

    common /grid/ dx, dy

    write(number,'(i4.4)') file
    file_name = "Rslt"//number
    open(unit=1, file=file_name//'.plt')

    write(1,*) 'TITLE="', time, '"'
    write(1,*) 'VARIABLES = "x","y","rho","vx","vy","Pre","alpha","grad"'
    write(1,*) "ZONE I=", NX, " J=", NY, " F=POINT"

    ! Compute total density over full domain including ghosts
    do j = -ghostp, NY+ghostp
        do i = -ghostp, NX+ghostp
            density(i,j) = alpha_rho_one(i,j) + alpha_rho_two(i,j)
        enddo
    enddo

    ! Density gradient components (central differences)
    do j = 1, NY
        do i = 1, NX
            density_gradx(i,j) = (density(i+1,j) - density(i-1,j))/(2.0d0*dx)
            density_grady(i,j) = (density(i,j+1) - density(i,j-1))/(2.0d0*dy)
        enddo
    enddo

    ! Schlieren-like gradient magnitude
    do j = 1, NY
        do i = 1, NX
            density_grad_mag(i,j) = dsqrt(density_gradx(i,j)**2.0d0 + density_grady(i,j)**2.0d0)
        enddo
    enddo

    dg_max = MAXVAL(density_grad_mag)

    do j = 1, NY
        do i = 1, NX
            write(1,'(9F25.8)') x(i), y(j), density(i,j), u_vel(i,j), v_vel(i,j),  &
                                 pressure(i,j), alpha_one(i,j), density_grad_mag(i,j), &
                                 dexp(density_grad_mag(i,j)/dg_max)
        enddo
    enddo

    close(1)

end subroutine output


! =============================================================================
! Initial conditions
!   Shock-droplet interaction test case:
!     - Left region (x < 2.89 cm): post-shock gas state
!     - Circular droplet centred at (4 cm, 3.7 cm), radius 1.1 cm: water
!     - Remaining domain: ambient air
!   EOS parameters: water (stiffened gas), air (ideal gas)
! =============================================================================

subroutine initialconditions(x, y, alpha_rho_one, alpha_rho_two, u_vel, v_vel,  &
                              pressure, alpha_one, NX, NY, test_case, t_end,     &
                              ghostp, n_eqn, gamma_one, gamma_two, pi_one, pi_two, dx, dy)

    implicit none

    integer          :: i, j
    integer          :: NX, NY, ghostp, n_eqn, test_case
    double precision :: x(-ghostp:NX+ghostp), y(-ghostp:NY+ghostp), dx, dy

    double precision :: t_end

    double precision :: density(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: pressure(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: u_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: v_vel(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: sound(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision :: alpha_rho_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_rho_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_one(-ghostp:NX+ghostp,-ghostp:NY+ghostp)
    double precision :: alpha_two(-ghostp:NX+ghostp,-ghostp:NY+ghostp)

    double precision :: gamma_one, gamma_two, pi_one, pi_two
    double precision :: gamma_big, pi_big
    double precision :: gamma, p_inf, dR, fsm

    double precision, parameter :: pi = acos(-1.0d0)

    ! End time
    t_end = 67d-6

    do i = 1, NX
        do j = 1, NY

            if (x(i) .lt. 2.89d0/100.0d0) then

                ! Post-shock region (left of incident shock)
                alpha_rho_one(i,j) = 1.0d0*(1.0d-8)
                alpha_rho_two(i,j) = 3.7579d0
                u_vel(i,j)         = 574.57d0
                v_vel(i,j)         = 0.0d0
                pressure(i,j)      = 6.6189d5
                alpha_one(i,j)     = 1.0d-8

            elseif ( ((x(i) - (4.0d0/100.0d0))**2.0d0                &
                   + (y(j) - (3.7d0/100.0d0))**2.0d0)                 &
                   .le. (1.1d0/100.0d0)**2.0d0 ) then

                ! Water droplet
                alpha_rho_one(i,j) = 1000.0d0
                alpha_rho_two(i,j) = 1.0d-8
                u_vel(i,j)         = 0.0d0
                v_vel(i,j)         = 0.0d0
                pressure(i,j)      = 1.013125d5
                alpha_one(i,j)     = (1.0d00 - 1.0d-8)

            else

                ! Ambient air
                alpha_rho_one(i,j) = 1.0d-8
                alpha_rho_two(i,j) = 1.17d0
                u_vel(i,j)         = 0.0d0
                v_vel(i,j)         = 0.0d0
                pressure(i,j)      = 1.013125d5
                alpha_one(i,j)     = 1.0d-8

            endif

            ! EOS constants: phase 1 = water (stiffened), phase 2 = air (ideal)
            gamma_one = 6.12d0
            gamma_two = 1.4d0
            pi_one    = 3.43d8
            pi_two    = 0.0d0

        enddo
    enddo

end subroutine initialconditions
