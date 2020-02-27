! 2020 Computer simulations of anisotropic structures in magnetorheological elastomers


subroutine write_start(n_part, d_cyl, d_part, h_cyl, phi)
    integer n_part
    real d_cyl, d_part, h_cyl, phi
    print '(a, i0, a)', 'Simulate ', n_part, ' particles in a cylinder:'
    print '(a, f0.2, a, f4.2, a)', 'd_cyl = ', d_cyl, ' = ', 1000.0*d_cyl*d_part, ' mm'
    print '(a, f0.2, a, f4.2, a)', 'h_cyl = ', h_cyl, ' = ', 1000.0*h_cyl*d_part, ' mm'
    print '(a, f4.1, a/)', 'at volume concentration of the particles: phi = ', phi*100.0, ' (v/v)%'
end subroutine write_start


! writing the table header to text files
subroutine write_vmd_title(n_part, h_cyl)
    integer n_part
    real h_cyl
    write(1, '(a, 3f10.6)') 'unitcell', h_cyl, h_cyl, h_cyl
    write(1, '(a, i0, a/)') 'atom 0:', n_part - 1, ' radius 0.5 name Z'
end subroutine write_vmd_title


! writting x, y, z coordinates of the particles to text files
subroutine write_vmd(n_part, x, y, z)
    integer i, n_part
    real x(n_part), y(n_part), z(n_part)
    write(1, '(a)') 'timestep ordered'
    do i = 1, n_part
        write(1, '(3f20.6)') x(i), y(i), z(i)
    end do
    write(1, '(a)') ''
end subroutine write_vmd
