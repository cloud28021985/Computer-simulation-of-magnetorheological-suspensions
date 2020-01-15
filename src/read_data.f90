! 2019 Computer simulations of anisotropic structures in magnetorheological elastomers


subroutine read_sigma(n_ext, chi, m_s, mu_0, pi, sigma_vs_gamma, sigma_vs_gamma_sat)
    integer i, j, n_ext
    real chi, m_s, mu_0, pi, res, res_sat, sigma_mean_kilo_si, sigma_mean_sat_kilo_si,  &
    sigma_mean(n_ext + 1), sigma_mean_sat(n_ext + 1), &
    gamma(n_ext + 1, 4), sigma(n_ext + 1, 4), sigma_sat(n_ext + 1, 4)
    real h(2) /65e3, 95e3/ ! magnetic field
    character(256) sigma_vs_gamma(4), sigma_vs_gamma_sat(4), sigma_mean_vs_gamma_sat

    character(256) sigma_mean_vs_gamma(2) /'data/sigma_mean_vs_gamma_h_1.txt', &
                                           'data/sigma_mean_vs_gamma_h_2.txt'/

    sigma_mean_vs_gamma_sat = 'data/sigma_mean_vs_gamma_sat.txt'

    print '(a/)', 'Reading sigma data'

    do j = 1, 4
        open(1, file = sigma_vs_gamma(j))
        open(2, file = sigma_vs_gamma_sat(j))
        do i = 1, n_ext + 1
            read (1, '(2f20.6)') gamma(i, j), sigma(i, j)
            read (2, '(2f20.6)') gamma(i, j), sigma_sat(i, j)
        end do
        close(1)
        close(2)
    end do

    do i = 1, n_ext + 1
        res = 0.0
        res_sat = 0.0
        do j = 1, 4
            res = res + sigma(i, j)
            res_sat = res_sat + sigma_sat(i, j)
        end do
        sigma_mean(i) = res / 4.0
        sigma_mean_sat(i) = res_sat / 4.0
    end do

    do j = 1, 2
        open(3, file = sigma_mean_vs_gamma(j))
        write(3, '(2a20)') '# gamma', 'sigma_kilo, Pa'
        do i = 1, n_ext + 1
            call worker_sigma_kilo_si(chi, h(j), mu_0, pi, sigma_mean(i), sigma_mean_kilo_si)
            write(3, '(2f20.6)') gamma(i, 1), sigma_mean_kilo_si
        end do
        close(3)
    end do

    open(4, file = sigma_mean_vs_gamma_sat)
    write(4, '(2a20)') '# gamma', 'sigma_sat_kilo, Pa'
    do i = 1, n_ext + 1
        call worker_sigma_sat_kilo_si(m_s, mu_0, pi, sigma_mean_sat(i), sigma_mean_sat_kilo_si)
        write(4, '(2f20.6)') gamma(i, 1), sigma_mean_sat_kilo_si
    end do
    close(4)
end subroutine read_sigma
















