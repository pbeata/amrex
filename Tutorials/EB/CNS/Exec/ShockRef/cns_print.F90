subroutine cns_printdata(level, time, lo, hi, u, ulo, uhi, dx, prob_lo, procNum, levNum, boxNum, stepNum) bind(C, name="cns_printdata")
  use amrex_fort_module, only : rt => amrex_real
  use cns_physics_module, only : gamma, cv
  use cns_module, only : center, nvar, urho, umx, umy, umz, ueden, ueint, utemp
  use probdata_module, only : p0,p1,rho0,rho1,v0,v1,x1

  implicit none

  integer, intent(in) :: level, lo(3), hi(3), ulo(3), uhi(3)
  real(rt), intent(in) :: time
  real(rt), intent(inout) :: u(ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3),nvar)
  real(rt), intent(in) :: dx(3), prob_lo(3)
  
  
  integer, intent(in) :: procNum, levNum, boxNum, stepNum
  integer :: i, j, k, count, ierr
  double precision :: x, y, z, xmax
  character (len=128) :: proc, lev, box, step, filename

  xmax = 3.5d0

  write(proc, '(I0)') procNum
  write(lev, '(I0)') levNum
  write(box, '(I0)') boxNum
  write(step, '(I0)') stepNum

  ! filename = 'error_' // trim(proc) // '.log'
  ! open(unit=1, file=filename, status='unknown')

  ! filename = 'output/proc' // trim(proc) // '_level' // trim(lev) // '_box' // trim(box) // '_step' // trim(step) // '.csv'
  filename = 'output/proc' // trim(proc) // '_level' // trim(lev) // '_box' // trim(box) // '.csv'

  open(unit=100, file=filename, status='unknown', iostat=ierr)
  ! if (ierr.ne.0) write(1, *) "FILE OPEN ERROR: ", ierr, filename

  count = 1

  do k = lo(3), hi(3)
    z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)

    do j = lo(2), hi(2)
      y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)

      do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

        if (x .gt. xmax) then
          write (100, '(I0, ", ", I0, ", ", I0, 4(",", ES24.16))') i, j, k, x, y, z, u(i,j,k,urho)
          count = count + 1
        end if

      end do
    end do
  end do

  ! close(1)
  close(100)

end subroutine cns_printdata
