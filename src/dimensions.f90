! 2020 Computer simulations of anisotropic structures in magnetorheological elastomers


subroutine worker_sigma_kilo_si(chi, h, mu_0, pi, sigma, sigma_kilo_si)
    real chi, h, mu_0, pi, sigma, sigma_kilo_si
    sigma_kilo_si = (3.0*pi*mu_0) / (16.0*1000.0) * pow_2(chi*h / (3.0 + chi)) * sigma
end subroutine worker_sigma_kilo_si


subroutine worker_sigma_sat_kilo_si(m_s, mu_0, pi, sigma, sigma_kilo_si)
    real m_s, mu_0, pi, sigma, sigma_kilo_si
    sigma_kilo_si = pi*mu_0*pow_2(m_s)*sigma / (48.0*1000.0)
end subroutine worker_sigma_sat_kilo_si


subroutine worker_t_micro_si(chi, eta, h, mu_0, t, t_micro_si)
    real chi, eta, h, mu_0, t, t_micro_si
    t_micro_si = 1000.0 * (16.0*eta) / mu_0 * pow_2((3.0 + chi) / (chi*h)) * t 
end subroutine worker_t_micro_si


subroutine worker_m_kilo_si(chi, h, m, m_kilo_si)
    real chi, h, m, m_kilo_si
    m_kilo_si = 3.0*chi*h*m / (1000.0*(3.0 + chi))
end subroutine worker_m_kilo_si
