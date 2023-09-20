program test
use iso_fortran_env, only: r8 => real64

use molecular_orbit
use read_fch
use advanced_fun

real(kind=r8), parameter :: ang2br = 1.8897260_r8 ! the unit in fch document is bohr 

integer, parameter :: iounit = 12

real(kind = r8), dimension(:,:), allocatable :: x_array
real(kind = r8), dimension(:,:), allocatable :: y_array
real(kind = r8), dimension(:,:), allocatable :: z_array 
real(kind = r8), dimension(:,:), allocatable :: rho, ELF_a, ELF_t

integer :: molecule_index

type(molecule) :: mol
 
character(len = 80) :: systemname,filename, outputfile, positionfile, numstring
real(kind = r8) :: x0, y0, z0 !the origin of the coordinate system
real(kind = r8) :: xmin, xmax
real(kind = r8) :: ymin, ymax
!real(kind = r8) :: zmin, zmax
integer :: resolution
integer :: iter, jter
!integer :: xiter, yiter, ziter

character(len = 80) :: written_format

!integrate on sphere
real(kind = r8) :: rmax
integer :: resolution_theta, resolution_phi, resolution_r
real(kind = r8), dimension(:,:), allocatable :: xyz_set
real(kind = r8), dimension(:), allocatable :: integration

procedure(adv_fun), pointer :: fun_ptr => null()




!-------------input parameters
systemname = "Au13opt6d10f"
molecule_index = 121
!*-------------paremeter for the span_xy
x0 = 0.0; y0 = 0.0; z0 = 0.0    !the origin of the coordinate system
xmin = -6.0 * ang2br; xmax = 6.0 * ang2br
ymin = -6.0 * ang2br; ymax = 6.0 * ang2br
resolution = 100


!*-------------integrate on sphere
resolution_theta = 40
resolution_phi = 80
resolution_r = 30
rmax = 6.0_r8 * ang2br

write(numstring, '(I0)') molecule_index
filename = trim(systemname) // ".fch"
outputfile = trim(systemname) // trim(numstring) // ".dat"
positionfile = trim(systemname) // trim(numstring) // "_position.dat"

!-------------input parameters

! the output format of the data
write(written_format, '(A,I0,A)')  '(', resolution, 'ES17.7E3)'

!!x_array and y_array have been allocated in span_xy
!!they can NOT be allocated twice!
!-------------allocation
allocate (x_array(resolution, resolution))
allocate (y_array(resolution, resolution))
allocate (z_array(resolution, resolution))
allocate (rho(resolution, resolution))
allocate (ELF_a(resolution, resolution))

!call span_xy(xmin, xmax, ymin, ymax, resolution, x_array, y_array)

!todo: build a module to generate a plane or 3d grid

!-------------read the data from the fch file
open(unit=iounit, file=filename, status='old', action='read')
call read_basis_func(iounit, mol, x0, y0, z0)
call read_molecular_orbit(iounit, mol)
close(iounit)

!todo: build a module to calculate the density for all kind of shell (R U NO RO)

!rho = 2*mol%den_a(x_array, y_array, z_array, molecule_index)
!ELF_a = ELF_alpha(mol, x_array, y_array, z_array)
!ELF_b = ELF_beta(mol, x_array, y_array, z_array)


!rho = 2*mol%den_a(x_array, y_array, z_array, molecule_index)

!ELF_t = ELF_total(mol, x_array, y_array, z_array, molecule_index)

!open(unit=iounit, file=outputfile, status='new', action='write')
!write(iounit, written_format) ((ELF_t(jter, iter), iter = 1, resolution), jter = 1, resolution)
!close(iounit)




!-------------integrate on sphere
allocate(integration(resolution_r))
fun_ptr => ELF_total

integration = sphere_integrate(rmax, & 
            resolution_r, resolution_theta, resolution_phi, &
            mol, molecule_index, fun_ptr)

open(unit=iounit, file="sphere_integrationg.txt", status='new', action='write')
write(iounit, '(1ES17.7E3)') (integration(iter), iter = 1, resolution_r)
close(iounit)

!-------------deallocation
deallocate(x_array)
deallocate(y_array)
deallocate(z_array)
deallocate(rho)
deallocate(ELF_a)
deallocate(integration)

write(* , *) "done!"


end program test







