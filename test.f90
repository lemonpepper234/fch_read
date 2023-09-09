module general_fun
    use iso_fortran_env, only: r8 => real64
    implicit none
    real(kind=r8) :: pi = 3.141592653589793_r8
   
    contains
   
       integer function factorial(n)
        integer, intent(in) :: n
        integer :: iter
   
        if (n < 0) stop "routine factorial: negative number input"
        factorial = 1
        do iter = 2, n
         factorial = factorial * iter
        end do
   
       end function factorial

       function linspace(a, b, n_elements)
        real(kind = r8), intent(in) :: a, b
        integer, intent(in), optional :: n_elements
        real(kind = r8), dimension(:), allocatable :: linspace

        real(kind = r8) :: dx
        integer :: i
        integer :: n
        integer :: ierr
        
        if (present(n_elements)) then
          if (n_elements <=1 ) then
              print*, "linspace procedure: Error: wrong value of n_elements, use an n_elements > 1"
              stop
          end if
          n=n_elements
          else
          n=100
        end if
    
        allocate(linspace(n), stat=ierr)
        if (ierr /= 0) then
          print*, "linspace procedure: Fatal Error, Allocation failed in linspace function"
          stop
        end if

        dx = (b - a) / (n - 1)
        linspace = (/ (i*dx + a, i = 1, n -1) /)

       end function linspace

       subroutine span_xy(xmin, xmax, ymin, ymax, resolution, x_array, y_array)
        real(kind = r8), intent(in) :: xmin, xmax
        real(kind = r8), intent(in) :: ymin, ymax
        real(kind = r8), dimension(:,:), allocatable :: x_array
        real(kind = r8), dimension(:,:), allocatable :: y_array
        integer, intent(in) :: resolution

        integer :: i

        do i = 1, resolution
          x_array(i, :) = linspace(xmin, xmax, resolution)
          y_array(:, i) = linspace(ymin, ymax, resolution)
        end do

       end subroutine span_xy
   
       subroutine locate_label(iounit, label)
        integer, intent(in) :: iounit
        character(len=*), intent(in) :: label
        character(len=80) :: char_tmp
        integer :: ierr
   
        rewind(iounit)
        do while (.true.)
         read(iounit, '(A80)', iostat=ierr) char_tmp
         if ( ierr /= 0) stop "routinte locate_label: label not found"
         if ( index(adjustl(char_tmp), label) == 1 ) exit
        end do
        backspace(iounit)
        !keep the position of the pointer for the document
       end subroutine locate_label
   
end module general_fun

module wave_fun
    use iso_fortran_env, only: r8 => real64
    use general_fun, f => factorial
    implicit none
    
     !  gtf1    gtf2     gtf3  |  gtf4      gtf5     gtf6  |  gtf7     gtf8     gtf9  |  gtf10
     !          basis_fun1(X)               basis_fun2(Y)              basis_fun3(Z)  |  basis_fun4(S)
     !                                      shell p                                      shell s
     !each basis_fun in the same shell has the same parameters
     !each gtf in the same contracted set has the same i j k (type) but different coefficients
        type :: gtf
         ! gtf = N * (x - x0)^i * (y - y0)^j * (z - z0)^k * exp( -zeta * r^2 )
         real(kind=r8) :: x0, y0, z0
         integer :: lx, ly, lz ! i, j, k 
         real(kind=r8) :: zeta
         real(kind=r8) :: Norm ! N: normalization factor
    
         contains
          procedure :: initialize => gtf_init
          procedure :: value => gtf_value
          procedure :: partialx => gtf_derivx
          procedure :: partialy => gtf_derivy
          procedure :: partialz => gtf_derivz
    
        end type gtf
    
        !the combination of gtfs (coefficients included)
        !one basis_fun contains multiple contracted gtfs
        type :: basis_fun
         type(gtf), dimension(:), allocatable :: gtfarray
         real(kind=r8), dimension(:), allocatable :: gtfcoeff
    
         contains
          procedure :: value => basis_value
          procedure :: partialx => basis_derivx
          procedure :: partialy => basis_derivy
          procedure :: partialz => basis_derivz
    
        end type basis_fun
    
     contains
    
    !-------------------------------------------------------function related to type gtf
        ! initialize the type gtf
        subroutine gtf_init( gtfType, locate, lx, ly, lz, zeta )
         class(gtf) :: gtfType !here gtf will change, so it is not intent(in)
         integer, intent(in) :: lx, ly, lz
         real(kind=r8), intent(in) :: zeta
         real(kind=r8), dimension(3), intent(in) :: locate
        
         gtfType%lx = lx
         gtfType%ly = ly
         gtfType%lz = lz
    
         gtfType%zeta = zeta
    
         gtfType%x0 = locate(1)
         gtfType%y0 = locate(2)
         gtfType%z0 = locate(3)
    
         gtfType%Norm = (2.*zeta/pi)**0.75 * sqrt( (8.*zeta)**(lx+ly+lz) &
                     * f(lx)*f(ly)*f(lz) / (f(2*lx)*f(2*ly)*f(2*lz)) )
       
        end subroutine gtf_init
    
        ! value of the gtf
        elemental real(kind = r8) function gtf_value(gtfType, x, y, z)
         class(gtf), intent(in) :: gtfType
         real(kind = r8), intent(in) :: x, y, z
    
         associate(dx => x - gtfType%x0, &
                   dy => y - gtfType%y0, &
                   dz => z - gtfType%z0)
         gtf_value = gtfType%Norm * dx**gtfType%lx &
                     * dy**gtfType%ly * dz**gtfType%lz &
                     * exp( -gtfType%zeta * (dx**2 + dy**2 + dz**2) )
         end associate
        
        end function gtf_value
    
        ! partial derivative alone x axis of the gtf
        elemental real(kind = r8) function gtf_derivx(gtfType, x, y, z)
         class(gtf), intent(in) :: gtfType
         real(kind = r8), intent(in) :: x, y, z
         
         associate(dx => x - gtfType%x0, &
                   dy => y - gtfType%y0, &
                   dz => z - gtfType%z0)
         
         if (gtfType%lx == 0) then
            gtf_derivx = -2.0*gtfType%zeta * dx &
                            * gtfType%Norm * dy**gtfType%ly * dz**gtfType%lz &
                            * exp( -gtfType%zeta * (dx**2 + dy**2 + dz**2) )
         else
            gtf_derivx = (gtfType%lx * dx**(gtfType%lx - 1) - 2.0 * gtfType%zeta * dx**(gtfType%lx + 1)) & 
                            * gtfType%Norm * dy**gtfType%ly * dz**gtfType%lz &
                            * exp( -gtfType%zeta * (dx**2 + dy**2 + dz**2) )
         end if
    
         end associate
    
        end function gtf_derivx
    
        ! partial derivative alone y axis of the gtf
        elemental real(kind = r8) function gtf_derivy(gtfType, x, y, z)
         class(gtf), intent(in) :: gtfType
         real(kind = r8), intent(in) :: x, y, z
    
         associate(dx => x - gtfType%x0, &
                   dy => y - gtfType%y0, &
                   dz => z - gtfType%z0)
         
         if (gtfType%ly == 0) then
            gtf_derivy = -2.0*gtfType%zeta * dy &
                            * gtfType%Norm * dx**gtfType%lx * dz**gtfType%lz &
                            * exp( -gtfType%zeta * (dx**2 + dy**2 + dz**2) )
         else
            gtf_derivy = (gtfType%ly * dy**(gtfType%ly - 1) - 2.0 * gtfType%zeta * dy**(gtfType%ly + 1)) & 
                            * gtfType%Norm * dx**gtfType%lx * dz**gtfType%lz &
                            * exp( -gtfType%zeta * (dx**2 + dy**2 + dz**2) )
         end if
    
         end associate
        
        end function gtf_derivy
    
        ! partial derivative alone z axis of the gtf
        elemental real(kind = r8) function gtf_derivz(gtfType, x, y, z)
         class(gtf), intent(in) :: gtfType
         real(kind = r8), intent(in) :: x, y, z
    
         associate(dx => x - gtfType%x0, &
                   dy => y - gtfType%y0, &
                   dz => z - gtfType%z0)   
         
         if (gtfType%lz == 0) then
            gtf_derivz = -2.0*gtfType%zeta * dz &
                            * gtfType%Norm * dx**gtfType%lx * dy**gtfType%ly &
                            * exp( -gtfType%zeta * (dx**2 + dy**2 + dz**2) )
         else
            gtf_derivz = (gtfType%lz * dz**(gtfType%lz - 1) - 2.0 * gtfType%zeta * dz**(gtfType%lz + 1)) & 
                            * gtfType%Norm * dx**gtfType%lx * dy**gtfType%ly &
                            * exp( -gtfType%zeta * (dx**2 + dy**2 + dz**2) )
         end if
    
         end associate
    
        end function gtf_derivz
    !-------------------------------------------------------function related to type gtf
    
    !-------------------------------------------------------function related to type basis_fun
        ! initialize the type basis_fun
        elemental real(kind = r8) function basis_value(basis_funType, x, y, z)
         class(basis_fun), intent(in) :: basis_funType
         real(kind = r8), intent(in) :: x, y, z
    
         basis_value = sum( basis_funType%gtfcoeff(:) * basis_funType%gtfarray(:)%value(x, y, z) )
    
        end function basis_value
    
        ! partial derivative alone x axis of the basis_fun
        elemental real(kind = r8) function basis_derivx(basis_funType, x, y, z)
          class(basis_fun), intent(in) :: basis_funType
          real(kind = r8), intent(in) :: x, y, z
    
          basis_derivx = sum( basis_funType%gtfcoeff(:) * basis_funType%gtfarray(:)%partialx(x, y, z) )
    
        end function basis_derivx
    
        ! partial derivative alone y axis of the basis_fun
        elemental real(kind = r8) function basis_derivy(basis_funType, x, y, z)
          class(basis_fun), intent(in) :: basis_funType
          real(kind = r8), intent(in) :: x, y, z
    
          basis_derivy = sum( basis_funType%gtfcoeff(:) * basis_funType%gtfarray(:)%partialy(x, y, z) )
    
        end function basis_derivy
    
        ! partial derivative alone z axis of the basis_fun
        elemental real(kind = r8) function basis_derivz(basis_funType, x, y, z)
          class(basis_fun), intent(in) :: basis_funType
          real(kind = r8), intent(in) :: x, y, z
    
          basis_derivz = sum( basis_funType%gtfcoeff(:) * basis_funType%gtfarray(:)%partialz(x, y, z) )
    
        end function basis_derivz
    
    !-------------------------------------------------------function related to type basis_fun
    
end module wave_fun

module molecular_orbit
    use iso_fortran_env, only: r8 => real64
    use wave_fun

  type :: molecule
    integer :: num_elec ! number of electrons
    integer :: num_alpha ! number of alpha electrons
    integer :: num_beta ! number of beta electrons
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

end module molecular_orbit

module read_fch
    use iso_fortran_env, only: r8 => real64
    use molecular_orbit
    implicit none
    
     !---------------------------------------------------------
     ! local variables
     ! but are global for submodule subroutines
     integer :: num_basis
     integer :: num_contra_shl
     integer :: num_primitive
     integer :: num_mo
     integer, dimension(:), allocatable :: shl_type
     integer, dimension(:), allocatable :: shl_ctr_odr
     real(kind=r8), dimension(:), allocatable :: prim_exp
     real(kind=r8), dimension(:), allocatable :: ctr_cff
     real(kind=r8), dimension(:), allocatable :: SP_ctr_cff
     real(kind=r8), dimension(:), allocatable :: shl_coord
     !---------------------------------------------------------
    
    contains
     subroutine read_basis_func(iounit, mol, x0, y0, z0)
      integer, intent(in) :: iounit
      type(molecule), intent(inout) :: mol
      real(kind=r8), intent(in) :: x0, y0, z0 !the origin of the coordinate system
    !-------------------------------------------------------
       integer, dimension(-3:3), parameter :: type2nbss = [7,5,4,1,3,6,10]
       ! reference Sobereva's codes: http://sobereva.com/55
       integer, dimension(-3:3,30) :: s2f
       integer, dimension(20), parameter :: type2lx = [0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1]
       integer, dimension(20), parameter :: type2ly = [0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1]
       integer, dimension(20), parameter :: type2lz = [0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1]
       !-------------------------------------------------------
       character(len=80) :: char_tmp
       integer :: lx, ly, lz
       real(kind=r8), dimension(3) :: locate
       integer :: ibss, ishl, ibsshl, ictr, iprim
       !-------------------------------------------------------
       s2f(0,1) = 1
       s2f(-1,1:4) = [1,2,3,4]
       s2f(1,1:3) = [2,3,4]
       s2f(2,1:6) = [5,6,7,8,9,10]
       s2f(3,1:10) = [11,12,13,17,14,15,18,19,16,20]
       !-------------------------------------------------------
       !num_contra_shl---(expand the conpoent of shell, such as X Y Z for p)----num_basis 
       !num_contra_shl---(expand the contraction part, 
       !while doesn't expand the component part 
       !cause the parameters of the components of the same type are the same)---num_primitive and ctr_cff
    
       
        ! read number of basis functions
        call locate_label(iounit, 'Number of basis functions')
        read(iounit, '(50X,I11)') num_basis
    
        ! read number of contracted shells
        call locate_label(iounit, 'Shell types')
        read(iounit, '(50X,I11)') num_contra_shl
       
        ! read type of each shell
        allocate( shl_type(num_contra_shl) )
        read(iounit, *) shl_type ! just utilize list-directed i/o
    
        ! read contraction order of each shell   
        read(iounit, *) ! skip one line
        allocate( shl_ctr_odr(num_contra_shl) )
        read(iounit, *) shl_ctr_odr
    
        ! read number of contracted shells
        call locate_label(iounit, 'Primitive exponents')
        read(iounit, '(50X,I11)') num_primitive
    
        ! read exponents of primitives (GTF)
        allocate( prim_exp(num_primitive) )
        read(iounit, *) prim_exp
    
        ! read contraction coefficients
        read(iounit, *) ! skip one line
        allocate( ctr_cff(num_primitive) )
        read(iounit, *) ctr_cff
    
        ! read SP(P) contraction coefficients if exist
        read(iounit, '(A80)') char_tmp
        if ( index(adjustl(char_tmp), &
             'P(S=P) Contraction coefficients') == 1) then
          allocate( SP_ctr_cff(num_primitive) )
          read(iounit, *) SP_ctr_cff
        end if
    
        ! read coordinates of contracted shells
        call locate_label(iounit, 'Coordinates of each shell ')
        read(iounit, *) 
        allocate( shl_coord(3*num_contra_shl) )
        read(iounit, *) shl_coord
    
        !-------------------------------------------------------
    
       allocate( mol%basis(num_basis) )
       ibss = 1; iprim = 1
    
       ! loop over shells
       do ishl = 1, num_contra_shl
    
          !loop over basis functions in the same shell
          do ibsshl = 1, type2nbss(shl_type(ishl))
          
          associate( basis => mol%basis(ibss), constrct_order => shl_ctr_odr(ishl) )
    
          allocate( basis%gtfarray(constrct_order) )
          allocate( basis%gtfcoeff(constrct_order) )
    
          locate = shl_coord(ishl*3-2:ishl*3) - (/x0, y0, z0/)
          lx = type2lx( s2f(shl_type(ishl), ibsshl) )
          ly = type2ly( s2f(shl_type(ishl), ibsshl) )
          lz = type2lz( s2f(shl_type(ishl), ibsshl) )
    
             ! loop over contracted GTFs in the same components of the same shell
             do ictr = 1, constrct_order
              if ( shl_type(ishl) == -1 .and. ibsshl /= 1 ) then
                ! sp type shell, and the last three GTFs are p-functions
                basis%gtfcoeff(ictr) = SP_ctr_cff(ictr+iprim-1) 
              else
                ! sp type shell, the first is s-function
                basis%gtfcoeff(ictr) = ctr_cff(ictr+iprim-1)
              end if
    
              !SP type shell, only subsitute the last three GTFs(p-functions)
    
              call basis%gtfarray(ictr)%initialize(locate, lx, ly, lz, prim_exp(ictr+iprim-1))
             end do
    
          end associate
          
          ibss = ibss + 1
          
          end do
    
          iprim = iprim + shl_ctr_odr(ishl)
          
       end do
    
     end subroutine read_basis_func
    
     subroutine read_molecular_orbit(iounit, mol)
       integer, intent(in) :: iounit
       type(molecule), intent(inout) :: mol
       character(len=80) :: char_tmp
       character(len=10) :: jobtype, method
       integer :: imo, ibss
    
       ! check wave function type
       rewind(iounit)
       read(iounit, *) char_tmp
       if (index(adjustl(char_tmp), 'isNO') == 1) then
          mol%wf_type = 'NO'
       else
          read(iounit, '(A10,A10)') jobtype, method
          ! confused by methods of which the name begins with letter 'O'
          if (index(trim(method), 'RO') == 1) then 
            mol%wf_type = 'RO'
          else if (index(trim(method), 'R') == 1) then
            mol%wf_type = 'R'
          else if (index(trim(method), 'U') == 1) then
            mol%wf_type = 'U'
          end if
        end if
    
       ! read the number of electrons (alpha/beta)
        call locate_label(iounit, 'Number of electrons')
        read(iounit, '(50X,I11)') mol%num_elec
        read(iounit, '(50X,I11)') mol%num_alpha
        read(iounit, '(50X,I11)') mol%num_beta
       
       ! read the number of MOs
        call locate_label(iounit, 'Alpha Orbital Energies')
        read(iounit, '(50X,I11)') num_mo
       
        ! if wave function is saved as natural orbitals
        ! 'Orbital Energies' is occupation numbers
        if ( trim(mol%wf_type) == 'NO' ) then
          allocate( mol%num_occ(num_mo) )
          read(iounit, *) mol%num_occ
        end if
    
        ! read MO coefficients
        allocate( mol%amo_coeff(num_mo, num_basis) )
        call locate_label(iounit, 'Alpha MO coefficients')
        read(iounit, *) ! skip one line
        read(iounit, *) ( ( mol%amo_coeff(imo,ibss), &
                         ibss=1,num_basis ), imo=1,num_mo )
    
        ! if is unrestricted wave function, read beta MO coefficients
        ! beta orbitals have no difference with alpha orbitals 
        ! in unrestricted natural orbital
        if ( mol%wf_type == 'U' ) then
          allocate( mol%bmo_coeff(num_mo, num_basis) )
          call locate_label(iounit, 'Beta MO coefficients')
          read(iounit, *) ! skip one line
          read(iounit, *) ( ( mol%bmo_coeff(imo,ibss), &
                           ibss=1,num_basis ), imo=1,num_mo )
        end if
    
        if ( .not. allocated(mol%bmo_coeff) ) then
          allocate( mol%bmo_coeff(num_mo, num_basis) )
          mol%bmo_coeff = mol%amo_coeff
        end if
       
    
     end subroutine read_molecular_orbit
    
    
end module read_fch


program test
use iso_fortran_env, only: r8 => real64
use read_fch

real(kind=r8), parameter :: ang2br = 1.8897260_r8 ! the unit in fch document is bohr 

integer, parameter :: iounit = 12
logical :: out_exist, pos_exist

real(kind = r8), dimension(:,:), allocatable :: x_array
real(kind = r8), dimension(:,:), allocatable :: y_array
real(kind = r8), dimension(:,:), allocatable :: z_array 
real(kind = r8), dimension(:,:), allocatable :: rho

type(molecule) :: mol

character(len = 80) :: filename, outputfile, positionfile
real(kind = r8) :: x0, y0, z0 !the origin of the coordinate system
real(kind = r8) :: xmin, xmax
real(kind = r8) :: ymin, ymax
real(kind = r8) :: zmin, zmax
integer :: resolution
integer :: iter, jter
integer :: xiter, yiter, ziter

character(len = 80) :: written_format

filename = "CH3_6-31Gdp.fch"
outputfile = "CH3_6-31Gdp.dat"
positionfile = "CH3_6-31Gdp_position.dat"
x0 = 0.0; y0 = 0.0; z0 = 0.0
xmin = -2.0 * ang2br; xmax = 2.0 * ang2br
ymin = -2.0 * ang2br; ymax = 2.0 * ang2br
resolution = 1000

! the output format of the data
write(written_format, '(A,I0,A)')  '(', resolution, 'ES15.7)'

!-------------allocation
allocate (x_array(resolution, resolution))
allocate (y_array(resolution, resolution))
allocate (z_array(resolution, resolution))
allocate (rho(resolution, resolution))

call span_xy(xmin, xmax, ymin, ymax, resolution, x_array, y_array)

z_array = 0.0_r8


open(unit=iounit, file=filename, status='old', action='read')
call read_basis_func(iounit, mol, x0, y0, z0)
call read_molecular_orbit(iounit, mol)
close(iounit)

rho = mol%den_a(x_array, y_array, z_array, 2) + mol%den_b(x_array, y_array, z_array, 2)



open(unit=iounit, file=outputfile, status='new', action='write')
write(iounit, written_format) ((rho(iter, jter), jter = 1, resolution),iter = 1, resolution)
close(iounit)





!write(*,'(30ES15.7)') mol%basis(:)%value(0.0_r8, 0.0_r8, 0.0_r8)
!write(*,'(30ES15.7)') mol%amo_coeff(2, :)
!write(*,'(30ES15.7)') mol%amo_coeff(2, :) * mol%basis(:)%value(0.0_r8,0.0_r8,0.0_r8)
!write(*,'(30ES15.7)') mol%amo_coeff(2, :)
!write(*,'(30ES15.7)') mol%bmo_coeff(2, :)

!write(*,*) mol%basis(1)%value(0.0_r8, 0.0_r8, 0.0_r8)
!write(*,*) mol%basis(1)%gtfcoeff(:)
!write(*,*) mol%basis(1)%gtfarray(:)%Norm
!write(*,*) mol%basis(1)%gtfarray(:)%zeta
!write(*,*) f(0)



!-------------deallocation
deallocate(x_array)
deallocate(y_array)
deallocate(z_array)
deallocate(rho)

end program test







