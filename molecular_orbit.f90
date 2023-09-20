module molecular_orbit
  use iso_fortran_env, only: r8 => real64
  use wave_fun

 type :: molecule
  integer :: num_elec ! number of electrons
  integer :: num_alpha ! number of alpha electrons
  integer :: num_beta ! number of beta electrons
  integer :: num_orbit
  character(len=2) :: wf_type ! wave function type: 'R', 'U', 'RO', 'NO'
  type(basis_fun), dimension(:), allocatable :: basis
  real(kind=r8), dimension(:), allocatable :: num_occ
  real(kind=r8), dimension(:,:), allocatable :: amo_coeff ! C(imo,ibss)
  real(kind=r8), dimension(:,:), allocatable :: bmo_coeff

   contains
    procedure :: den_a => density_alpha
    procedure :: den_b => density_beta
    procedure :: den_n => density_natural
    procedure :: den_partialx_a => density_partialx_alpha
    procedure :: den_partialy_a => density_partialy_alpha
    procedure :: den_partialz_a => density_partialz_alpha
    procedure :: den_partialx_b => density_partialx_beta
    procedure :: den_partialy_b => density_partialy_beta
    procedure :: den_partialz_b => density_partialz_beta
    procedure :: den_partialx_n => density_partialx_natural
    procedure :: den_partialy_n => density_partialy_natural
    procedure :: den_partialz_n => density_partialz_natural
    procedure :: den_lapal_a => density_lapal_alpha
    procedure :: den_lapal_b => density_lapal_beta
 end type

 contains
 
 !-------------------------------------------functions to calculate the density
 !the occupied number of alpha orbit euqals to 1
 elemental real(kind = r8) function density_alpha(moleculeType, x, y, z, num_molecule)! string can not be the input paramenter of the elemental function
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den

  temp_den = sum(moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  density_alpha = temp_den**2 

 end function density_alpha

 !the occupied number of beta orbit euqals to 1
 elemental real(kind = r8) function density_beta(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den

  temp_den = sum(moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  density_beta = temp_den**2 

 end function density_beta

 !the occupied number of naturala orbit is not integer--->moleculeType%num_occ(num_molecule)
 elemental real(kind = r8) function density_natural(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
 
  temp_den = sum(moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  density_natural = temp_den**2 * moleculeType%num_occ(num_molecule)

 end function density_natural
 !-------------------------------------------functions to calculate the density

 !-------------------------------------------functions to calculate the partial derivative of the density of alpha orbit
 elemental real(kind = r8) function density_partialx_alpha(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
  real(kind = r8) :: temp_den_partialx

  temp_den = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  temp_den_partialx = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%partialx(x, y, z) )

  density_partialx_alpha = 2.0 * temp_den * temp_den_partialx

 end function density_partialx_alpha

 elemental real(kind = r8) function density_partialy_alpha(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
  real(kind = r8) :: temp_den_partialy

  temp_den = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  temp_den_partialy = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%partialy(x, y, z) )

  density_partialy_alpha = 2.0 * temp_den * temp_den_partialy

 end function density_partialy_alpha

 elemental real(kind = r8) function density_partialz_alpha(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
  real(kind = r8) :: temp_den_partialz

  temp_den = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  temp_den_partialz = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%partialz(x, y, z) )

  density_partialz_alpha = 2.0 * temp_den * temp_den_partialz

 end function density_partialz_alpha
 !-------------------------------------------functions to calculate the partial derivative of the density of alpha orbit

 !-------------------------------------------functions to calculate the partial derivative of the density of beta orbit
 elemental real(kind = r8) function density_partialx_beta(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
  real(kind = r8) :: temp_den_partialx

  temp_den = sum( moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  temp_den_partialx = sum( moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%partialx(x, y, z) )

  density_partialx_beta = 2.0 * temp_den * temp_den_partialx

 end function density_partialx_beta

 elemental real(kind = r8) function density_partialy_beta(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
  real(kind = r8) :: temp_den_partialy

  temp_den = sum( moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  temp_den_partialy = sum( moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%partialy(x, y, z) )

  density_partialy_beta = 2.0 * temp_den * temp_den_partialy

 end function density_partialy_beta

 elemental real(kind = r8) function density_partialz_beta(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
  real(kind = r8) :: temp_den_partialz

  temp_den = sum( moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  temp_den_partialz = sum( moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%partialz(x, y, z) )

  density_partialz_beta = 2.0 * temp_den * temp_den_partialz

 end function density_partialz_beta
 !-------------------------------------------functions to calculate the partial derivative of the density of beta orbit

 !-------------------------------------------functions to calculate the partial derivative of the density of natural orbit
 elemental real(kind = r8) function density_partialx_natural(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
  real(kind = r8) :: temp_den_partialx

  num_basis = size(moleculeType%basis)
  temp_den = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  temp_den_partialx = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%partialx(x, y, z) )

  density_partialx_natural = 2.0 * temp_den * temp_den_partialx * moleculeType%num_occ(num_molecule)

 end function density_partialx_natural

 elemental real(kind = r8) function density_partialy_natural(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
  real(kind = r8) :: temp_den_partialy

  num_basis = size(moleculeType%basis)
  temp_den = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  temp_den_partialy = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%partialy(x, y, z) )

  density_partialy_natural = 2.0 * temp_den * temp_den_partialy * moleculeType%num_occ(num_molecule)

 end function density_partialy_natural

 elemental real(kind = r8) function density_partialz_natural(moleculeType, x, y, z, num_molecule)
  class(molecule), intent(in) :: moleculeType
  real(kind = r8), intent(in) :: x, y, z
  integer, intent(in) :: num_molecule
  real(kind = r8) :: temp_den
  real(kind = r8) :: temp_den_partialz

  num_basis = size(moleculeType%basis)
  temp_den = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%value(x, y, z) )
  temp_den_partialz = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%partialz(x, y, z) )

  density_partialz_natural = 2.0 * temp_den * temp_den_partialz * moleculeType%num_occ(num_molecule)

 end function density_partialz_natural

 !-------------------------------------------functions to calculate the partial derivative of the density of natural orbit

 !-------------------------------------------functions to calculate the lapalcian function of orbit function

 elemental real(kind = r8) function density_lapal_alpha(moleculeType, x, y, z, num_molecule)
 class(molecule), intent(in) :: moleculeType
 real(kind = r8), intent(in) :: x, y, z
 integer, intent(in) :: num_molecule
 real(kind = r8) :: temp_den_partialx, temp_den_partialy, temp_den_partialz

 temp_den_partialx = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%partialx(x, y, z) )
 temp_den_partialy = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%partialy(x, y, z) )
 temp_den_partialz = sum( moleculeType%amo_coeff(num_molecule, :) * moleculeType%basis(:)%partialz(x, y, z) )

 density_lapal_alpha = temp_den_partialx**2 + temp_den_partialy**2 + temp_den_partialz**2

 end function density_lapal_alpha

 elemental real(kind = r8) function density_lapal_beta(moleculeType, x, y, z, num_molecule)
 class(molecule), intent(in) :: moleculeType
 real(kind = r8), intent(in) :: x, y, z
 integer, intent(in) :: num_molecule
 real(kind = r8) :: temp_den_partialx, temp_den_partialy, temp_den_partialz

 temp_den_partialx = sum( moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%partialx(x, y, z) )
 temp_den_partialy = sum( moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%partialy(x, y, z) )
 temp_den_partialz = sum( moleculeType%bmo_coeff(num_molecule, :) * moleculeType%basis(:)%partialz(x, y, z) )

 density_lapal_beta = temp_den_partialx**2 + temp_den_partialy**2 + temp_den_partialz**2

 end function density_lapal_beta

 !-------------------------------------------functions to calculate the lapalcian function of orbit function



end module molecular_orbit