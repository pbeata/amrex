subroutine writeout(lo, hi, u, ulo, uhi, dx, prob_lo, procNum, levNum, boxNum) bind(C, name="writeout")

  ! mandatory declarations for calling Fortran from AMReX
  use amrex_fort_module, only : amrex_real
  integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3)
  real(amrex_real), intent(inout) :: u( ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3) )

  ! my declarations
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer, intent(in) :: procNum, levNum, boxNum
  integer :: i, j, k, count
  double precision :: x, y, z
  character (len=128) :: proc, lev, box, filename

  ! main work kernel
  ! (let's assume 2D case only for now)
  count = 1

  write(proc, '(I0)') procNum
  write(lev, '(I0)') levNum
  write(box, '(I0)') boxNum
  
  !filename = 'output/proc' // trim(proc) // '_box' // trim(box) // '.txt'
  filename = 'output/proc' // trim(proc) // '_lev' // trim(lev) // '_box' // trim(box) // '.csv'
  open(unit=100, file=filename, status='unknown')
  

  !do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
      do i=lo(1),hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

        ! a better output format for CSV file
        write (100, '(I0, ", ", I0, 3(",", ES24.16))') i, j, x, y, u(i,j,0)

        count = count + 1
      end do
    end do
  !end do

  close(100)

end subroutine writeout