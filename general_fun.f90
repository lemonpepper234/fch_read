module general_fun
  use iso_fortran_env, only: r8 => real64, i8 => int64
  implicit none
  real(kind=r8) :: pi = 3.141592653589793_r8

  contains

  ! warning: possible to overflow
     integer(kind = i8) function factorial(n)
      integer, intent(in) :: n
      integer :: iter
 
      if (n < 0) then
       print *, "routine factorial: negative number input, n = ", n ; 
       stop "Terminating due to negative input to factorial."
      end if

      factorial = 1_i8
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
      linspace = (/ (i*dx + a, i = 0, n -1) /)

     end function linspace
     
     !!x_array and y_array have already been allocated
     subroutine span_xy(xmin, xmax, ymin, ymax, resolution, x_array, y_array)
      real(kind = r8), intent(in) :: xmin, xmax
      real(kind = r8), intent(in) :: ymin, ymax
      integer, intent(in) :: resolution

      integer :: i

      real(kind = r8), dimension(:,:), allocatable, intent(out) :: x_array
      real(kind = r8), dimension(:,:), allocatable, intent(out) :: y_array

      allocate(x_array(resolution, resolution))
      allocate(y_array(resolution, resolution))

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