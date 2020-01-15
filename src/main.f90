! 2019 Computer simulations of anisotropic structures in magnetorheological elastomers


subroutine main_alrorythm(n_part, n_ext, n_int, b, d_cyl, d_part, eps_st, h_cyl, phi, pi, r_cut, time_step, &
    time_step_forcecap, vtf, cross_section, sigma_vs_gamma, model, rank)
    include 'mpif.h'
    integer i, model, n_ext, n_int, n_part, n_up, n_down, n_free, rank, &
    bottom_wall(n_part), upper_wall(n_part), free_part(n_part)
    real b, d_cyl, d_part, delta_gamma, delta_t, eps_st, finish_time, gamma, gamma_max, h_cyl, phi, pi, r_cut, &
    r_min, sigma, start_time, t, time_step, time_step_forcecap, &
    x(n_part),   y(n_part),   z(n_part), &
    m_x(n_part), m_y(n_part), m_z(n_part), &
    f_x(n_part), f_z(n_part)
    character(256) cross_section, vtf, sigma_vs_gamma

    t = 0.0

    if(rank == 0) then
        call write_start(n_part, d_cyl, d_part, h_cyl, phi)
    end if

    ! generate $n_part particles at random positions
    call worker_random(n_part, d_cyl, h_cyl, x, y, z)

    ! open VTF file
    open(1, file = vtf)
    call write_vmd_title(n_part, h_cyl)

    ! WARMUP INTEGRATION
    if(rank == 0) then
        print '(a)', 'Warmup:'
    end if

    call cpu_time(start_time)

    ! compute the minimal distance between two particles
    call mindist(n_part, x, y, z, r_min)

    i = 0
    do while(r_min < 1.0)
        call integrate_forcecap(n_part, d_cyl, h_cyl, time_step_forcecap, x, y, z)

        ! warmup criterion
        call mindist(n_part, x, y, z, r_min)

        i = i + 1
    end do

    call cpu_time(finish_time)

    if(rank == 0) then
        print '(a, i0, a, f0.6/)', 'step = ', i, ', act_min_dist = ', r_min
    end if

    delta_t = finish_time - start_time

    if(rank == 0) then
        print '(a, f4.2, a/)', 'warmup runtime = ', delta_t, ' sec'
    end if

    call worker_m_initial(n_part, m_x, m_y, m_z)

    ! MAIN INTEGRATION
    if(rank == 0) then
        print '(a)', 'Main integration:'
    end if

    call cpu_time(start_time)

    do i = 0, n_ext
        if(rank == 0) then
            print '(1a1, a, f5.1, a, $)', char(13), 'percent complete ', 100.0*i / n_ext, '%'
        end if

        if(model == 1) then
            call worker_m(n_part, b, x, y, z, m_x, m_y, m_z)
        end if

        ! output configuration
        call write_vmd(n_part, x, y, z)

        call integrate(n_part, n_int, d_cyl, h_cyl, eps_st, r_cut, time_step, x, y, z, m_x, m_y, m_z)
        t = t + time_step*n_int
    end do

    call cpu_time(finish_time)
    delta_t = finish_time - start_time

    if(rank == 0) then
        print '(1a1, a, f0.2, a, $)', char(13), 'main runtime = ', delta_t, ' sec = '
        print '(f0.2, a, f0.2, a/)', delta_t / 60.0, ' min = ', delta_t / 3600.0, ' hours'
    end if

    ! open cross section file
    open(4, file = cross_section)
    call write_cross_section_title(n_part, h_cyl)

    call write_cross_section(n_part, x, y, z)
    close(4)

    ! SHEAR INTEGRATION
    if(rank == 0) then
        print '(a)', 'Shear integration:'
    end if

    ! open shear file
    open(3, file = sigma_vs_gamma)

    call cpu_time(start_time)
    call detecting_wall_part(n_part, n_up, n_down, n_free, bottom_wall, upper_wall, free_part, h_cyl, y)

    gamma = 0.0
    gamma_max = 0.5
    delta_gamma = gamma_max / (1.0*n_ext*n_int)

    do i = 0, n_ext
        if(rank == 0) then
            print '(1a1, a, f5.1, a, $)', char(13), 'percent complete ', 100.0*i / n_ext, '%'
        end if

        if(model == 1) then
            call worker_m(n_part, b, x, y, z, m_x, m_y, m_z)
        end if

        call worker_f_shear(n_part, n_up, n_free, upper_wall, free_part, eps_st, r_cut, x, y, z, m_x, m_y, m_z, &
        f_x, f_z)

        call worker_sigma(n_part, n_up, upper_wall, d_cyl, pi, sigma, x, z, f_x, f_z)

        ! output configuration
        call write_vmd(n_part, x, y, z)

        write(3, '(2f20.6)') gamma, sigma

        call integrate_shear(n_part, n_up, n_int, upper_wall, n_free, free_part, d_cyl, delta_gamma, eps_st, &
        h_cyl, r_cut, time_step, x, y, z, m_x, m_y, m_z)

        gamma = gamma + delta_gamma*n_int
    end do

    call cpu_time(finish_time)
    delta_t = finish_time - start_time

    if(rank == 0) then
        print '(1a1, a, f0.2, a, $)', char(13), 'shear runtime = ', delta_t, ' sec = '
        print '(f0.2, a, f0.2, a/)', delta_t / 60.0, ' min = ', delta_t / 3600.0, ' hours'
    end if

    close(1)
    close(3)
end subroutine main_alrorythm


subroutine plan_simulate(n_part, n_ext, n_int, b, d_cyl, d_part, eps_st, h_cyl, phi, pi, r_cut, time_step, &
    time_step_forcecap, sigma_vs_gamma, sigma_vs_gamma_sat, rank, size)
    include 'mpif.h'
    integer i, n_part, n_ext, n_int, rank, size
    real b, d_cyl, d_part, eps_st, h_cyl, phi, pi, r_cut, time_step, time_step_forcecap
    character(256) sigma_vs_gamma(4), sigma_vs_gamma_sat(4)

    character(256) vtf(4) /'vmd/r_vs_t_experiment_1.vtf', &
                           'vmd/r_vs_t_experiment_2.vtf', &
                           'vmd/r_vs_t_experiment_3.vtf', &
                           'vmd/r_vs_t_experiment_4.vtf'/

    character(256) vtf_sat(4) /'vmd/r_vs_t_experiment_sat_1.vtf', &
                               'vmd/r_vs_t_experiment_sat_2.vtf', &
                               'vmd/r_vs_t_experiment_sat_3.vtf', &
                               'vmd/r_vs_t_experiment_sat_4.vtf'/

    character(256) cross_section(4) /'vmd/cross_section_1.vtf', &
                                     'vmd/cross_section_2.vtf', &
                                     'vmd/cross_section_3.vtf', &
                                     'vmd/cross_section_4.vtf'/

    character(256) cross_section_sat(4) /'vmd/cross_section_sat_1.vtf', &
                                         'vmd/cross_section_sat_2.vtf', &
                                         'vmd/cross_section_sat_3.vtf', &
                                         'vmd/cross_section_sat_4.vtf'/

    if(rank == 0) then
        print '(a/)', 'Main alrorythm for various phi'
    end if

    do i = 1, size
        if(rank == i - 1) then
            call main_alrorythm(n_part, n_ext, n_int, b, d_cyl, d_part, eps_st, h_cyl, phi, pi, r_cut, time_step, &
            time_step_forcecap, vtf(i), cross_section(i), sigma_vs_gamma(i), 1, rank)

            call main_alrorythm(n_part, n_ext, n_int, b, d_cyl, d_part, eps_st, h_cyl, phi, pi, r_cut, time_step, &
            time_step_forcecap, vtf_sat(i), cross_section_sat(i), sigma_vs_gamma_sat(i), 0, rank)
        end if
    end do
end subroutine plan_simulate


subroutine plan_read(n_ext, chi, delta_t, m_s, mu_0, pi, sigma_vs_gamma, sigma_vs_gamma_sat, rank)
    integer n_ext, rank
    real chi, m_s, mu_0, pi
    double precision delta_t
    character(256) sigma_vs_gamma(4), sigma_vs_gamma_sat(4)
    if(rank == 0) then
        print '(a//)', '----------------------------'
        print '(a, f0.2, a, f0.2, a, $)', 'Total runtime = ', delta_t, ' sec = ', delta_t / 60.0, ' min = '
        print '(f0.2, a/)', delta_t / 3600.0, ' hours'
        call read_sigma(n_ext, chi, m_s, mu_0, pi, sigma_vs_gamma, sigma_vs_gamma_sat)
        print '(a/)', 'Finished' ! terminate program
    end if
end subroutine plan_read


program main
    include 'mpif.h'
    integer, parameter :: &

    !n_part = 2, &
    n_part = 40, &
    !n_part = 1000, & ! number of the particles in the system

    n_ext = 40 ! number of external nodes
    integer ierr, n_int, rank, size
    real, parameter :: &
    chi = 1000.0, & ! magnetic susceptibility of the particle
    d_part = 5e-6, & ! diameter of the particle
    eps_st = 2.0, & ! parameter of sterical interaction
    m_s = 1655e3, & ! suturate magnetization of the particles
    pi = 4.0*atan(1.0), &
    phi = 0.3, & ! volume concentration of the particles
    time_step = 0.001, & ! set up the integrator time step
    time_step_forcecap = 0.01
    real b, mu_0, r_cut
    double precision delta_t, finish_time, start_time

    character(256) sigma_vs_gamma(4) /'data/sigma_vs_gamma_experiment_1.txt', &
                                      'data/sigma_vs_gamma_experiment_2.txt', &
                                      'data/sigma_vs_gamma_experiment_3.txt', &
                                      'data/sigma_vs_gamma_experiment_4.txt'/

    character(256) sigma_vs_gamma_sat(4) /'data/sigma_vs_gamma_experiment_sat_1.txt', &
                                          'data/sigma_vs_gamma_experiment_sat_2.txt', &
                                          'data/sigma_vs_gamma_experiment_sat_3.txt', &
                                          'data/sigma_vs_gamma_experiment_sat_4.txt'/

    form_cyl = 1.0
    !form_cyl = 5.7 ! form-factor of cylinder

    n_int = 2000
    !n_int = 22000 ! number of internal nodes

    b = chi / (8.0*(3.0 + chi))
    mu_0 = 4*pi*1e-7 ! vacuum permeability
    r_cut = 2.0 ** (1.0 / 6.0)

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, size, ierr)
    call mpi_comm_rank(mpi_comm_world, rank, ierr)

    if(n_part == 2) then
        n_int = 60
        form_cyl = 1.0
    end if

    if(rank == 0) then
        print '(a//)', '----------------------------'
    end if

    start_time = mpi_wtime(ierr)

    call worker_cylinder(n_part, d_cyl, form_cyl, h_cyl, phi)

    call plan_simulate(n_part, n_ext, n_int, b, d_cyl, d_part, eps_st, h_cyl, phi, pi, r_cut, time_step, &
    time_step_forcecap, sigma_vs_gamma, sigma_vs_gamma_sat, rank, size)

    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)
    finish_time = mpi_wtime(ierr)
    delta_t = finish_time - start_time

    call plan_read(n_ext, chi, delta_t, m_s, mu_0, pi, sigma_vs_gamma, sigma_vs_gamma_sat, rank)
end program main

! Total runtime = 20223.44 sec = 337.06 min = 5.62 hours

















