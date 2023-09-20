module advanced_fun
        use iso_fortran_env, only: r8 => real64, i8 => int64
        use general_fun
        use molecular_orbit
        implicit none

        !!elemental function can NOT be a procedure pointer!!
        abstract interface  
                real(kind = r8) function adv_fun(moleculeType, x, y, z, num_molecule)
                import :: r8
                import :: molecule
                class(molecule), intent(in) :: moleculeType
                real(kind = r8), intent(in) :: x, y, z
                integer, intent(in) :: num_molecule
                end function adv_fun
        end interface

        contains
 !------------------------the function to generate the points distributed on the sphere
        !!xyz_set has already been allocated!!
        !*xyz_set -> [resolution_theta*resolution_phi,5] the last two column are theta and phi
        subroutine span_sphere_xyz(r, resolution_theta, resolution_phi, xyz_set)
                real(kind = r8), intent(in) :: r
                integer, intent(in) :: resolution_theta, resolution_phi

                real(kind = r8), dimension(:,:), allocatable, intent(out) :: xyz_set

                real(kind = r8), dimension(:), allocatable :: theta_array, phi_array
                integer :: theta_iter, phi_iter

                allocate(theta_array(resolution_theta))
                allocate(phi_array(resolution_phi))

                allocate(xyz_set(resolution_theta * resolution_phi, 5))

                theta_array = linspace(0.0_r8, pi, resolution_theta)
                phi_array = linspace(0.0_r8, 2.0_r8 * pi, resolution_phi)

                do theta_iter = 1, resolution_theta
                        do phi_iter = 1, resolution_phi
                                xyz_set((theta_iter-1)*resolution_phi + phi_iter, :) = &
                                 (/ r * sin(theta_array(theta_iter)) * cos(phi_array(phi_iter)) , &
                                 r * sin(theta_array(theta_iter)) * sin(phi_array(phi_iter)) , &
                                 r * cos(theta_array(theta_iter)) , &
                                 theta_array(theta_iter), phi_array(phi_iter) /)
                        end do
                end do

        end subroutine span_sphere_xyz

        function sphere_integrate(rmax, resolution_r, resolution_theta, resolution_phi, &
                 moleculeType, num_molecule, fun_ptr) result(integration)
                real(kind = r8), intent(in) :: rmax
                integer, intent(in) :: resolution_r, resolution_theta, resolution_phi
                class(molecule), intent(in) :: moleculeType
                integer, intent(in) :: num_molecule

                procedure(adv_fun), pointer, intent(in) :: fun_ptr

                real(kind = r8), dimension(:), allocatable :: integration
                
                real(kind = r8), dimension(:), allocatable :: r_array
                real(kind = r8), dimension(:,:), allocatable :: xyz_set
                integer :: r_iter, point_iter

                allocate(r_array(resolution_r))
                allocate(integration(resolution_r))

                r_array = linspace(0.0_r8, rmax, resolution_r)
                integration = 0.0_r8

                do r_iter = 1, resolution_r
                        do point_iter = 1, resolution_theta * resolution_phi
                                call span_sphere_xyz(r_array(r_iter), & 
                                resolution_theta, resolution_phi, xyz_set)

                                integration(r_iter) = integration(r_iter) + &
                                fun_ptr(moleculeType, xyz_set(point_iter, 1), &
                                xyz_set(point_iter, 2), xyz_set(point_iter, 3), &
                                num_molecule) * &
                                r_array(r_iter)**2 * sin(xyz_set(point_iter, 4)) * &
                                (rmax / resolution_r) * &
                                ( pi / resolution_theta) * &
                                (2.0_r8 * pi / resolution_phi)

                        end do   
                end do

        end function sphere_integrate 

 !------------------------the function to calculate the electron localization function for alpha electron 
        elemental real(kind = r8) function ELF_alpha_plane(moleculeType, x, y, z, num_molecule)
                class(molecule), intent(in) :: moleculeType
                real(kind = r8), intent(in) :: x, y, z
                integer, intent(in) :: num_molecule
                
                real(kind = r8) :: rho 
                real(kind = r8) :: rho_partialx 
                real(kind = r8) :: rho_partialy 
                real(kind = r8) :: rho_partialz 
                real(kind = r8) :: kinetic_energy 
                real(kind = r8) :: fermion_kinetic_energy 
                real(kind = r8) :: boson_kinetic_energy
                real(kind = r8) :: uniform_kinetic_energy

  
                associate(kin => kinetic_energy, &
                        f_kin => fermion_kinetic_energy, &
                        b_kin => boson_kinetic_energy, &
                        u_kin => uniform_kinetic_energy)
                
                        rho = moleculeType%den_a(x, y, z, num_molecule)
                        rho_partialx = moleculeType%den_partialx_a(x, y, z, num_molecule)
                        rho_partialy = moleculeType%den_partialy_a(x, y, z, num_molecule)
                        rho_partialz = moleculeType%den_partialz_a(x, y, z, num_molecule)
                        kin = 0.5_r8 * moleculeType%den_lapal_a(x, y, z, num_molecule)
                
                        b_kin = 0.25_r8 * & 
                                (rho_partialx**2 + rho_partialy**2 + rho_partialz**2) / &
                                rho
                        
                        u_kin = 0.6_r8 * (6.0_r8 * pi**2)**(2.0_r8/3.0_r8) * rho**(5.0_r8/3.0_r8)
                        
                        f_kin = kin - b_kin + 1E-5_r8
                        
                        ELF_alpha_plane = 1.0_r8 / (1.0_r8 + (f_kin / u_kin)**2.0_r8)
                
                end associate
  
        end function ELF_alpha_plane

        elemental real(kind = r8) function ELF_beta_plane(moleculeType, x, y, z, num_molecule)
                class(molecule), intent(in) :: moleculeType
                real(kind = r8), intent(in) :: x, y, z
                integer, intent(in) :: num_molecule

                real(kind = r8) :: rho
                real(kind = r8) :: rho_partialx
                real(kind = r8) :: rho_partialy
                real(kind = r8) :: rho_partialz
                real(kind = r8) :: kinetic_energy
                real(kind = r8) :: fermion_kinetic_energy
                real(kind = r8) :: boson_kinetic_energy
                real(kind = r8) :: uniform_kinetic_energy

                associate(kin => kinetic_energy, &
                        f_kin => fermion_kinetic_energy, &
                        b_kin => boson_kinetic_energy, &
                        u_kin => uniform_kinetic_energy)

                        rho = moleculeType%den_b(x, y, z, num_molecule)
                        rho_partialx = moleculeType%den_partialx_b(x, y, z, num_molecule)
                        rho_partialy = moleculeType%den_partialy_b(x, y, z, num_molecule)
                        rho_partialz = moleculeType%den_partialz_b(x, y, z, num_molecule)
                        kin = 0.5_r8 * moleculeType%den_lapal_b(x, y, z, num_molecule)

                        b_kin = 0.25_r8 * &
                                (rho_partialx**2 + rho_partialy**2 + rho_partialz**2) / &
                                rho

                        u_kin = 0.6_r8 * (6.0_r8 * pi**2)**(2.0_r8/3.0_r8) * rho**(5.0_r8/3.0_r8)

                        f_kin = kin - b_kin + 1E-5_r8

                        ELF_beta_plane = 1.0_r8 / (1.0_r8 + (f_kin / u_kin)**2.0_r8)

                end associate

        end function ELF_beta_plane
  
        elemental real(kind = r8) function ELF_total_plane(moleculeType, x, y, z, num_molecule)
                class(molecule), intent(in) :: moleculeType
                real(kind = r8), intent(in) :: x, y, z
                integer, intent(in) :: num_molecule
        
                real(kind = r8) :: rho 
                real(kind = r8) :: rho_partialx 
                real(kind = r8) :: rho_partialy 
                real(kind = r8) :: rho_partialz 
                real(kind = r8) :: kinetic_energy 
                real(kind = r8) :: fermion_kinetic_energy 
                real(kind = r8) :: boson_kinetic_energy
                real(kind = r8) :: uniform_kinetic_energy
        
                associate(kin => kinetic_energy, &
                f_kin => fermion_kinetic_energy, &
                b_kin => boson_kinetic_energy, &
                u_kin => uniform_kinetic_energy)
        
        
                        rho =  moleculeType%den_a(x, y, z, num_molecule) + &
                                        moleculeType%den_b(x, y, z, num_molecule)
                        rho_partialx =  moleculeType%den_partialx_a(x, y, z, num_molecule) + &
                                        moleculeType%den_partialx_b(x, y, z, num_molecule)
                        rho_partialy =  moleculeType%den_partialy_a(x, y, z, num_molecule) + &
                                        moleculeType%den_partialy_b(x, y, z, num_molecule)
                        rho_partialz =  moleculeType%den_partialz_a(x, y, z, num_molecule) + &
                                        moleculeType%den_partialz_b(x, y, z, num_molecule)
                        kin =  0.5_r8 * moleculeType%den_lapal_a(x, y, z, num_molecule) + &
                                        0.5_r8 * moleculeType%den_lapal_b(x, y, z, num_molecule)

                        b_kin = 0.125_r8 * & 
                                (rho_partialx**2 + rho_partialy**2 + rho_partialz**2) / &
                                rho

                        u_kin = 0.3_r8 * (6.0_r8 * pi**2)**(2.0_r8/3.0_r8) * rho**(5.0_r8/3.0_r8)

                        f_kin = kin - b_kin + 1E-5_r8

                        ELF_total_plane = 1.0_r8 / (1.0_r8 + (f_kin / u_kin)**2.0_r8)
        
                end associate
  
        end function ELF_total_plane

 !------------------------the ELF suit for procedure pointer
        real(kind = r8) function ELF_total(moleculeType, x, y, z, num_molecule)
        class(molecule), intent(in) :: moleculeType
        real(kind = r8), intent(in) :: x, y, z
        integer, intent(in) :: num_molecule

        real(kind = r8) :: rho 
        real(kind = r8) :: rho_partialx 
        real(kind = r8) :: rho_partialy 
        real(kind = r8) :: rho_partialz 
        real(kind = r8) :: kinetic_energy 
        real(kind = r8) :: fermion_kinetic_energy 
        real(kind = r8) :: boson_kinetic_energy
        real(kind = r8) :: uniform_kinetic_energy

        associate(kin => kinetic_energy, &
        f_kin => fermion_kinetic_energy, &
        b_kin => boson_kinetic_energy, &
        u_kin => uniform_kinetic_energy)


                rho =  moleculeType%den_a(x, y, z, num_molecule) + &
                                moleculeType%den_b(x, y, z, num_molecule)
                rho_partialx =  moleculeType%den_partialx_a(x, y, z, num_molecule) + &
                                moleculeType%den_partialx_b(x, y, z, num_molecule)
                rho_partialy =  moleculeType%den_partialy_a(x, y, z, num_molecule) + &
                                moleculeType%den_partialy_b(x, y, z, num_molecule)
                rho_partialz =  moleculeType%den_partialz_a(x, y, z, num_molecule) + &
                                moleculeType%den_partialz_b(x, y, z, num_molecule)
                kin =  0.5_r8 * moleculeType%den_lapal_a(x, y, z, num_molecule) + &
                                0.5_r8 * moleculeType%den_lapal_b(x, y, z, num_molecule)

                b_kin = 0.125_r8 * & 
                        (rho_partialx**2 + rho_partialy**2 + rho_partialz**2) / &
                        rho

                u_kin = 0.3_r8 * (6.0_r8 * pi**2)**(2.0_r8/3.0_r8) * rho**(5.0_r8/3.0_r8)

                f_kin = kin - b_kin + 1E-5_r8

                ELF_total = 1.0_r8 / (1.0_r8 + (f_kin / u_kin)**2.0_r8)

        end associate

end function ELF_total

end module advanced_fun