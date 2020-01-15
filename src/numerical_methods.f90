! 2019 Computer simulations of anisotropic structures in magnetorheological elastomers


! squaring
function pow_2(x)
    real pow_2, x 
    pow_2 = x*x
end function pow_2


function pow_3(x)
    real pow_3, x 
    pow_3 = x*x*x
end function pow_3


function pow_4(x)
    real pow_4, x, x_2
    x_2 = x*x
    pow_4 = x_2*x_2
end function pow_4


function pow_5(x)
    real pow_5, x, x_2
    x_2 = x*x
    pow_5 = x_2*x_2*x
end function pow_5


function pow_6(x)
    real pow_6, x, x_2
    x_2 = x*x
    pow_6 = x_2*x_2*x_2
end function pow_6


subroutine integrate_forcecap(n_part, d_cyl, h_cyl, time_step_forcecap, x, y, z)
    integer i, n_part
    real d_cyl, h_cyl, time_step_forcecap, &
    x(n_part),   y(n_part),   z(n_part), &
    f_x(n_part), f_y(n_part), f_z(n_part)
    call worker_wall_forcecap(n_part, d_cyl, h_cyl, x, y, z, f_x, f_y, f_z)
    call forcecap(n_part, x, y, z, f_x, f_y, f_z)
    do i = 1, n_part
        x(i) = x(i) + time_step_forcecap*f_x(i)
        y(i) = y(i) + time_step_forcecap*f_y(i)
        z(i) = z(i) + time_step_forcecap*f_z(i)
    end do
end subroutine integrate_forcecap


subroutine integrate(n_part, n_int, d_cyl, h_cyl, eps_st, r_cut, time_step, x, y, z, m_x, m_y, m_z)
    integer i, j, n_int, n_part
    real d_cyl, h_cyl, eps_st, r_cut, time_step, &
    x(n_part),   y(n_part),   z(n_part), &
    f_x(n_part), f_y(n_part), f_z(n_part), &
    m_x(n_part), m_y(n_part), m_z(n_part)
    do i = 1, n_int
        call worker_wall(n_part, d_cyl, h_cyl, r_cut, x, y, z, f_x, f_y, f_z)
        call worker_f(n_part, r_cut, eps_st, x, y, z, m_x, m_y, m_z, f_x, f_y, f_z)
        do j = 1, n_part
            x(j) = x(j) + time_step*f_x(j)
            y(j) = y(j) + time_step*f_y(j)
            z(j) = z(j) + time_step*f_z(j)
        end do
    end do
end subroutine integrate


subroutine integrate_shear(n_part, n_up, n_int_shear, upper_wall, n_free, free_part, d_cyl, delta_gamma, eps_st, h_cyl, r_cut, &
time_step, x, y, z, m_x, m_y, m_z)
    integer n_part, n_up, n_int_shear, n_free, &
    free_part(n_part)
    real d_cyl, delta_gamma, eps_st, h_cyl, r_cut, time_step, &
    x(n_part),   y(n_part),   z(n_part), &
    f_x(n_part), f_y(n_part), f_z(n_part), &
    m_x(n_part), m_y(n_part), m_z(n_part), &
    upper_wall(n_part)
    do i = 1, n_int_shear
        call shear_upper_part(n_part, n_up, upper_wall, delta_gamma, d_cyl, h_cyl, x, z)
        call worker_wall(n_part, d_cyl, h_cyl, r_cut, x, y, z, f_x, f_y, f_z)
        call worker_f(n_part, r_cut, eps_st, x, y, z, m_x, m_y, m_z, f_x, f_y, f_z)
        do j = 1, n_free
            x(free_part(j)) = x(free_part(j)) + time_step*f_x(free_part(j))
            y(free_part(j)) = y(free_part(j)) + time_step*f_y(free_part(j))
            z(free_part(j)) = z(free_part(j)) + time_step*f_z(free_part(j))
        end do
    end do
end subroutine integrate_shear









