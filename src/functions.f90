! 2019 Computer simulations of anisotropic structures in magnetorheological elastomers


! define the system size
subroutine worker_cylinder(n_part, d_cyl, form_cyl, h_cyl, phi)
    integer n_part
    real d_cyl, form_cyl, h_cyl, phi
    d_cyl = (2.0*n_part*form_cyl / (3.0*phi)) ** (1.0/3.0)
    h_cyl = (2.0*n_part / (3.0*phi*pow_2(form_cyl))) ** (1.0/3.0)
end subroutine worker_cylinder


! random function with uniform distribution
subroutine worker_random(n_part, d_cyl, h_cyl, x, y, z)
    integer n_part
    real a, b, c, d_cyl, h_cyl, random_x, random_y, random_z, &
    x(n_part), y(n_part), z(n_part)
    a = (1.0 - d_cyl) / 2.0
    b = d_cyl - 1.0
    c = pow_2(d_cyl - 1.0) / 4.0
    i = 1
    do while (i <= n_part)
        call random_number(random_x)
        call random_number(random_y)
        call random_number(random_z)
        x(i) = a + b*random_x
        y(i) = 0.5 + (h_cyl - 1.0)*random_y
        z(i) = a + b*random_z
        if(pow_2(x(i)) + pow_2(z(i)) < c) then
            i = i + 1
        end if
    end do
end subroutine worker_random


subroutine worker_delta_x(x_i, x_j, y_i, y_j, z_i, z_j, delta_x, delta_y, delta_z)
    real x_i, x_j, y_i, y_j, z_i, z_j, delta_x, delta_y, delta_z
    delta_x = x_i - x_j
    delta_y = y_i - y_j
    delta_z = z_i - z_j
end subroutine worker_delta_x


subroutine mindist(n_part, x, y, z, r_min)
    integer i, j, n_part
    real delta_x, delta_y, delta_z, r_2, r_2_min, r_min, &
    x(n_part), y(n_part), z(n_part)
    call worker_delta_x(x(2), x(1), y(2), y(1), z(2), z(1), delta_x, delta_y, delta_z)
    r_2_min = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
    do i = 3, n_part
        do j = 1, i - 1
            call worker_delta_x(x(i), x(j), y(i), y(j), z(i), z(j), delta_x, delta_y, delta_z)
            r_2 = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
            if (r_2_min > r_2) then
                r_2_min = r_2
            end if
        end do
    end do
    r_min = sqrt(r_2_min)
end subroutine mindist


subroutine worker_wall_forcecap(n_part, d_cyl, h_cyl, x, y, z, f_x, f_y, f_z)
    integer i, n_part
    real a, b, d_cyl, h_cyl, r, r_2, &
    x(n_part),   y(n_part),   z(n_part), &
    f_x(n_part), f_y(n_part), f_z(n_part)
    a = h_cyl - 0.5
    b = pow_2(d_cyl - 1.0) / 4.0
    do i = 1, n_part
        if(y(i) < 0.5) then
            f_y(i) = 4.0
        else if(y(i) > a) then
            f_y(i) = - 4.0
        else
            f_y(i) = 0.0
        end if
        r_2 = pow_2(x(i)) + pow_2(z(i))
        if(r_2 > b) then
            r = sqrt(r_2)
            f_x(i) = - 4.0*x(i) / r
            f_z(i) = - 4.0*z(i) / r
        else
            f_x(i) = 0.0
            f_z(i) = 0.0
        end if
    end do
end subroutine worker_wall_forcecap


subroutine worker_wall(n_part, d_cyl, h_cyl, r_cut, x, y, z, f_x, f_y, f_z)
    integer i, n_part
    real a, b, c, d, d_cyl, h_cyl, r, r_2, r_cut, r_cyl, r_cyl_6, r_cyl_7, y_cyl, y_cyl_6, y_cyl_7, &
    x(n_part),   y(n_part),   z(n_part), &
    f_x(n_part), f_y(n_part), f_z(n_part)
    a = r_cut - 0.5
    b = h_cyl - r_cut + 0.5
    c = pow_2(0.5*d_cyl - r_cut + 0.5)
    do i = 1, n_part
        if(y(i) < a) then
            y_cyl = y(i) + 0.5
            y_cyl_6 = pow_6(y_cyl)
            y_cyl_7 = y_cyl*y_cyl_6
            f_y(i) = 2.0/ y_cyl_7 * (2.0 / y_cyl_6 - 1.0)
        else if(y(i) > b) then
            y_cyl = h_cyl + 0.5 - y(i)
            y_cyl_6 = pow_6(y_cyl)
            y_cyl_7 = y_cyl_6*y_cyl
            f_y(i) = 2.0/ y_cyl_7 * (1.0 - 2.0 / y_cyl_6)
        else
            f_y(i) = 0.0
        end if
        r_2 = pow_2(x(i)) + pow_2(z(i))
        if(r_2 > c) then
            r = sqrt(r_2)
            r_cyl = 0.5*d_cyl + 0.5 - r
            r_cyl_6 = pow_6(r_cyl)
            r_cyl_7 = r_cyl*r_cyl_6
            d = 2.0/ (r*r_cyl_7) * (1.0 - 2.0 / r_cyl_6)
            f_x(i) = d*x(i)
            f_z(i) = d*z(i)
        else
            f_x(i) = 0.0
            f_z(i) = 0.0
        end if
    end do
end subroutine worker_wall


subroutine forcecap(n_part, x, y, z, f_x, f_y, f_z)
    integer i, j, n_part
    real delta_x, delta_y, delta_z, force_x, force_y, force_z, r, r_2, &
    x(n_part),   y(n_part),   z(n_part), &
    f_x(n_part), f_y(n_part), f_z(n_part)
    do i = 1, n_part - 1
        do j = i + 1, n_part
            call worker_delta_x(x(i), x(j), y(i), y(j), z(i), z(j), delta_x, delta_y, delta_z)
            r_2 = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
            r = sqrt(r_2)
            if(r < 1.0) then
                force_x = delta_x / r
                force_y = delta_y / r
                force_z = delta_z / r
            else
                force_x = 0.0
                force_y = 0.0
                force_z = 0.0
            end if
            f_x(i) = f_x(i) + force_x
            f_y(i) = f_y(i) + force_y
            f_z(i) = f_z(i) + force_z
            f_x(j) = f_x(j) - force_x
            f_y(j) = f_y(j) - force_y
            f_z(j) = f_z(j) - force_z
        end do
    end do
end subroutine forcecap


! initial conditions for magnetization of the particles
subroutine worker_m_initial(n_part, m_x, m_y, m_z)
    integer i, n_part
    real m_x(n_part), m_y(n_part), m_z(n_part)
    do i = 1, n_part
        m_x(i) = 0.0
        m_y(i) = 1.0
        m_z(i) = 0.0
    end do
end subroutine worker_m_initial


subroutine worker_a(b, delta_x, delta_y, delta_z, a_xx, a_xy, a_xz, a_yy, a_yz, a_zz)
    real a_xx, a_xy, a_xz, a_yy, a_yz, a_zz, b, delta_x, delta_x_2, delta_y, delta_y_2, delta_z, delta_z_2, r, r_2, r_5
    delta_x_2 = delta_x*delta_x
    delta_y_2 = delta_y*delta_y
    delta_z_2 = delta_z*delta_z
    r_2 = delta_x_2 + delta_y_2 + delta_z_2
    r = sqrt(r_2)
    r_5 = r*r_2*r_2
    a_xx = b*(delta_y_2 + delta_z_2 - 2.0*delta_x_2) / r_5
    a_xy = - 3.0*b*delta_x*delta_y / r_5
    a_xz = - 3.0*b*delta_x*delta_z / r_5
    a_yy = b*(delta_x_2 + delta_z_2 - 2.0*delta_y_2) / r_5
    a_yz = - 3.0*b*delta_y*delta_z / r_5
    a_zz = b*(delta_x_2 + delta_y_2 - 2.0*delta_z_2) / r_5
end subroutine


! magnetization of the particles
subroutine worker_m(n_part, b, x, y, z, m_x, m_y, m_z)
    integer i, j, n_part
    real a_xx, a_xy, a_xz, a_yy, a_yz, a_zz, b, delta_x, delta_y, delta_z, error, norm, var_x, var_y, var_z, &
    x(n_part),   y(n_part),   z(n_part), &
    m_x(n_part), m_y(n_part), m_z(n_part), &
    p_x(n_part), p_y(n_part), p_z(n_part)
    error = 1e-5 ! error
    do
        do i = 1, n_part
            p_x(i) = m_x(i)
            p_y(i) = m_y(i)
            p_z(i) = m_z(i)
        end do
        do i = 1, n_part
            var_x = 0.0
            var_y = 0.0
            var_z = 0.0
            do j = 1, i - 1
                call worker_delta_x(x(i), x(j), y(i), y(j), z(i), z(j), delta_x, delta_y, delta_z)
                call worker_a(b, delta_x, delta_y, delta_z, a_xx, a_xy, a_xz, a_yy, a_yz, a_zz)
                var_y = var_y + a_yy*m_y(j) + a_xy*p_x(j) + a_yz*p_z(j)
                var_x = var_x + a_xy*m_y(j) + a_xx*m_x(j) + a_xz*p_z(j)
                var_z = var_z + a_yz*m_y(j) + a_xz*m_x(j) + a_zz*m_z(j)
            end do
            do j = i + 1, n_part
                call worker_delta_x(x(i), x(j), y(i), y(j), z(i), z(j), delta_x, delta_y, delta_z)
                call worker_a(b, delta_x, delta_y, delta_z, a_xx, a_xy, a_xz, a_yy, a_yz, a_zz)
                var_y = var_y + a_yy*p_y(j) + a_xy*p_x(j) + a_yz*p_z(j)
                var_x = var_x + a_xy*m_y(j) + a_xx*p_x(j) + a_xz*p_z(j)
                var_z = var_z + a_yz*m_y(j) + a_xz*m_x(j) + a_zz*p_z(j)
            end do
            m_x(i) = - var_x
            m_y(i) = 1.0 - var_y
            m_z(i) = - var_z
        end do
        norm = 0.0
        do i = 1, n_part
            norm = norm + pow_2(m_x(i) - p_x(i)) + pow_2(m_y(i) - p_y(i)) + pow_2(m_z(i) - p_z(i))
        end do
        if(sqrt(norm) < error) exit
    end do
end subroutine worker_m


! magnetic and sterical forces of interaction between particles
subroutine worker_f(n_part, r_cut, eps_st, x, y, z, m_x, m_y, m_z, f_x, f_y, f_z)
    integer i, j, n_part
    real a, b, c, delta_x, delta_y, delta_z, dot_m_m, dot_m_i_r, dot_m_j_r, eps_st, force_x, force_y, force_z, &
    force_m_x, force_m_y, force_m_z, force_wca_x, force_wca_y, force_wca_z, r, r_2, r_4, r_5, r_6, r_8, r_cut, &
    x(n_part),   y(n_part),   z(n_part), &
    m_x(n_part), m_y(n_part), m_z(n_part), &
    f_x(n_part), f_y(n_part), f_z(n_part)
    do i = 1, n_part - 1
        do j = i + 1, n_part
            call worker_delta_x(x(i), x(j), y(i), y(j), z(i), z(j), delta_x, delta_y, delta_z)
            r_2 = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
            r = sqrt(r_2)
            r_4 = r_2*r_2
            r_5 = r*r_4
            dot_m_m = m_x(i)*m_x(j) + m_y(i)*m_y(j) + m_z(i)*m_z(j)
            dot_m_i_r = m_x(i)*delta_x + m_y(i)*delta_y + m_z(i)*delta_z
            dot_m_j_r = m_x(j)*delta_x + m_y(j)*delta_y + m_z(j)*delta_z
            a = dot_m_i_r*dot_m_j_r / r_2
            b = dot_m_m - 5.0*a
            force_m_x = (delta_x*b + m_x(i)*dot_m_j_r + m_x(j)*dot_m_i_r) / r_5
            force_m_y = (delta_y*b + m_y(i)*dot_m_j_r + m_y(j)*dot_m_i_r) / r_5
            force_m_z = (delta_z*b + m_z(i)*dot_m_j_r + m_z(j)*dot_m_i_r) / r_5
            if(r > r_cut) then
                f_x(i) = f_x(i) + force_m_x
                f_y(i) = f_y(i) + force_m_y
                f_z(i) = f_z(i) + force_m_z
                f_x(j) = f_x(j) - force_m_x
                f_y(j) = f_y(j) - force_m_y
                f_z(j) = f_z(j) - force_m_z
            else
                r_6 = r_2*r_4
                r_8 = r_4*r_4
                c = 2.0*eps_st*(2.0 / r_6 - 1.0) / r_8
                force_wca_x = delta_x*c
                force_wca_y = delta_y*c
                force_wca_z = delta_z*c
                force_x = force_m_x + force_wca_x
                force_y = force_m_y + force_wca_y
                force_z = force_m_z + force_wca_z
                f_x(i) = f_x(i) + force_x
                f_y(i) = f_y(i) + force_y
                f_z(i) = f_z(i) + force_z
                f_x(j) = f_x(j) - force_x
                f_y(j) = f_y(j) - force_y
                f_z(j) = f_z(j) - force_z
            end if
        end do
    end do
end subroutine worker_f


subroutine worker_f_shear(n_part, n_up, n_free, upper_wall, free_part, eps_st, r_cut, x, y, z, m_x, m_y, m_z, f_x, f_z)
    integer i, j, n_part, n_up, n_free, &
    upper_wall(n_part), free_part(n_part)
    real a, b, c, delta_x, delta_y, delta_z, dot_m_m, dot_m_i_r, dot_m_j_r, eps_st, force_x, force_z, &
    force_m_x, force_m_z, force_wca_x, force_wca_z, r, r_2, r_4, r_5, r_6, r_8, r_cut, &
    x(n_part),   y(n_part),   z(n_part), &
    m_x(n_part), m_y(n_part), m_z(n_part), &
    f_x(n_part), f_z(n_part)
    f_x = 0.0
    f_z = 0.0
    do i = 1, n_up
        do j = 1, n_free
            call worker_delta_x(x(upper_wall(i)), x(free_part(j)), y(upper_wall(i)), y(free_part(j)), z(upper_wall(i)), &
            z(free_part(j)), delta_x, delta_y, delta_z)
            r_2 = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
            r = sqrt(r_2)
            r_4 = r_2*r_2
            r_5 = r*r_4
            r_6 = r_2*r_4
            r_8 = r_4*r_4
            dot_m_m = m_x(upper_wall(i))*m_x(free_part(j)) + m_y(upper_wall(i))*m_y(free_part(j)) + &
            m_z(upper_wall(i))*m_z(free_part(j))
            dot_m_i_r = m_x(upper_wall(i))*delta_x + m_y(upper_wall(i))*delta_y + m_z(upper_wall(i))*delta_z
            dot_m_j_r = m_x(free_part(j))*delta_x + m_y(free_part(j))*delta_y + m_z(free_part(j))*delta_z
            a = dot_m_i_r*dot_m_j_r / r_2
            b = dot_m_m - 5.0*a
            c = 2.0*eps_st*(2.0 / r_6 - 1.0) / r_8
            force_m_x = (delta_x*b + m_x(upper_wall(i))*dot_m_j_r + m_x(free_part(j))*dot_m_i_r) / r_5
            force_m_z = (delta_z*b + m_z(upper_wall(i))*dot_m_j_r + m_z(free_part(j))*dot_m_i_r) / r_5
            if(r > r_cut) then
                force_wca_x = 0.0
                force_wca_z = 0.0
            else
                force_wca_x = delta_x*c
                force_wca_z = delta_z*c
            end if
            force_x = force_m_x + force_wca_x
            force_z = force_m_z + force_wca_z
            f_x(upper_wall(i)) = f_x(upper_wall(i)) + force_x
            f_z(upper_wall(i)) = f_z(upper_wall(i)) + force_z
        end do
    end do
end subroutine worker_f_shear


subroutine detecting_wall_part(n_part, n_up, n_down, n_free, bottom_wall, upper_wall, free_part, h_cyl, y)
    integer i, n_part, n_up, n_down, n_free, &
    bottom_wall(n_part), upper_wall(n_part), free_part(n_part)
    real h_cyl, &
    y(n_part)
    n_down = 0
    n_up = 0
    n_free = 0
    bottom_wall = 0
    upper_wall = 0
    free_part = 0
    do i = 1, n_part
        if(y(i) < 1.0) then
            n_down = n_down + 1
            bottom_wall(n_down) = i
        else if(y(i) > h_cyl - 1.0) then
            n_up = n_up + 1
            upper_wall(n_up) = i
        else
            n_free = n_free + 1
            free_part(n_free) = i
        end if
    end do
    if(n_down + n_up + n_free .ne. n_part) then
        print '(a/, a/)', 'n_down + n_up + n_free .ne. n_part', 'stop'
        stop
    end if
end subroutine detecting_wall_part


subroutine shear_upper_part(n_part, n_up, upper_wall, delta_gamma, d_cyl, h_cyl, x, z)
    integer i, n_part, n_up, &
    upper_wall(n_part)
    real d_cyl, delta_gamma, h_cyl, &
    x(n_part), z(n_part)
    do i = 1, n_up
        z(upper_wall(i)) = z(upper_wall(i)) + x(upper_wall(i))*delta_gamma*h_cyl / d_cyl
        x(upper_wall(i)) = x(upper_wall(i)) - z(upper_wall(i))*delta_gamma*h_cyl / d_cyl
    end do
end subroutine shear_upper_part


subroutine worker_sigma(n_part, n_up, upper_wall, d_cyl, pi, sigma, x, z, f_x, f_z)
    integer i, j, n_part, n_up, &
    upper_wall(n_part)
    real d_cyl, delta, my_sum, pi, sigma, &
    x(n_part),   z(n_part), &
    f_x(n_part), f_z(n_part)
    delta = 1.0
    my_sum = 0.0
    do i = 1, n_up
        j = upper_wall(i)
        r = sqrt(pow_2(x(j)) + pow_2(z(j)))
        if(r > d_cyl / 2.0 - delta) then
            my_sum = my_sum + (z(j)*f_x(j) - x(j)*f_z(j)) / r
        end if
    end do
    sigma = 2.0*my_sum / (pi*d_cyl*delta)
end subroutine worker_sigma








