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