subroutine writeout(lo, hi, u, ulo, uhi, dx, prob_lo, boxNum) bind(C, name="writeout")

  ! mandatory declarations for calling Fortran from AMReX
  use amrex_fort_module, only : amrex_real
  integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3)
  real(amrex_real), intent(inout) :: u( ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3) )

  ! my declarations
  double precision, intent(in) :: dx(3), prob_lo(3)
  integer, intent(in) :: boxNum
  integer :: i, j, k, count
  double precision :: x, y, z
  character (len=64) :: filename

  ! main work kernel
  ! (let's assume 2D case only for now)
  count = 1

  if (boxNum .lt. 10) then
    write(filename, '( "output/test", I1, ".txt" )') boxNum
  else if (boxNum .lt. 100) then
    write(filename, '( "output/test", I2, ".txt" )') boxNum
  else if (boxNum .lt. 1000) then
    write(filename, '( "output/test", I3, ".txt" )') boxNum
  else
    ! error
  end if
  
  open(unit=100, file=filename, status='unknown')

  !do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
      do i=lo(1),hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
        write(100,*) i, j, x, y, u(i,j,0)
        count = count + 1
      end do
    end do
  !end do

  close(100)

end subroutine writeout