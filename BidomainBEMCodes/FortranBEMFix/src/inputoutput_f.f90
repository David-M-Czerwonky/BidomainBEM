module inputoutput
use globalconstants
implicit none
contains
!!!!!!! c suffix for complex data r suffix for real data 

subroutine initmeshfromfile(filename,trimesh,siz)
  implicit none
  integer(kind=spint),intent(in) :: siz
  character*(*),intent(in) :: filename
  type(trianglemesh),intent(out) :: trimesh
  real(kind=dp),dimension(siz) :: garb
  real(kind=dp),dimension(3) :: v1,v2
  integer(kind=spint) :: i,j
  call rreadfromdat(filename,siz,garb)
!!!!!allocate mesh arrays
  allocate(trimesh%p(3,trimesh%np))
i=1
j=3*trimesh%np
trimesh%p=reshape(dble(garb(i:j)),(/3,trimesh%np/))
  allocate(trimesh%t2p(3,trimesh%nt))
i=j+1
j=j+3*trimesh%nt
trimesh%t2p=reshape(idnint(garb(i:j)),(/3,trimesh%nt/))
  allocate(trimesh%reg(trimesh%nt))
i=j+1
j=j+trimesh%nt
trimesh%reg=garb(i:j)
  allocate(trimesh%nhat(3,trimesh%nt))
  allocate(trimesh%vol(trimesh%nt))
do i=1,trimesh%nt
v1=trimesh%p(:,trimesh%t2p(1,i))-trimesh%p(:,trimesh%t2p(3,i))
v2=trimesh%p(:,trimesh%t2p(2,i))-trimesh%p(:,trimesh%t2p(3,i))
trimesh%nhat(1,i)=v1(2)*v2(3)-v1(3)*v2(2)
trimesh%nhat(2,i)=v1(3)*v2(1)-v1(1)*v2(3)
trimesh%nhat(3,i)=v1(1)*v2(2)-v1(2)*v2(1)
trimesh%vol(i)=(trimesh%nhat(1,i)**2.0_dp+trimesh%nhat(2,i)**2.0_dp+trimesh%nhat(3,i)**2.0_dp)**0.5_dp
trimesh%nhat(:,i)=trimesh%nhat(:,i)/trimesh%vol(i)
trimesh%vol(i)=trimesh%vol(i)/2.0_dp;
end do

end subroutine initmeshfromfile

subroutine initcoilfromfile(filename,coil,siz)
  implicit none
  integer(kind=spint),intent(in) :: siz
  character*(*),intent(in) :: filename
  type(coilarr),intent(out) :: coil
  real(kind=dp),dimension(6*siz) :: garb
  allocate(coil%rs(3,siz)) 
  allocate(coil%js(3,siz))
  call rreadfromdat(filename,6*siz,garb)
  coil%rs=reshape(garb(1:3*siz),(/3,siz/))
  coil%js=reshape(garb(3*siz+1:6*siz),(/3,siz/))
end subroutine initcoilfromfile



!!!!!!! read file routines

 subroutine creadfromdat(filename,nx,exportdata)
    implicit none
    character*(*),intent(in) :: filename
    integer(kind=spint),intent(in) :: nx
    complex(kind=dp),dimension(nx),intent(out) :: exportdata
    real(kind=dp),dimension(nx,2) :: tempnum
    open(unit=11,file=filename,status='unknown',READONLY)
    read(11,'(E24.16E4)') tempnum
    write(*,*) tempnum(1:4,1)
    close(11)
    exportdata=cmplx(tempnum(:,1),tempnum(:,2),dp)
end subroutine creadfromdat	


subroutine rreadfromdat(filename,nx,exportdata)
  implicit none
  character*(*),intent(in) :: filename
  integer(kind=spint),intent(in) :: nx
  real(kind=dp),dimension(nx),intent(out) :: exportdata
  open(unit=11,file=filename,status='unknown',READONLY)
  read(11,'(E24.16E4)') exportdata
  write(*,*) exportdata(1:4)
  close(11)
end subroutine rreadfromdat




!!!!!!! write to file routines
subroutine cwritetofile(filename,array,s1)
    character*(*),intent(in) :: filename
    integer,intent(in) :: s1
    complex(kind=dp),dimension(:) :: array
   integer(kind=spint) :: i

   open(unit=10,file=filename,status='unknown') 
           do i=1,s1
             write(10,*) real(array(i)),imag(array(i))
          end do
    close(10)
end subroutine cwritetofile

  subroutine rwritetofile(filename,array,s1)
    character*(*),intent(in) :: filename
    integer,intent(in) :: s1
    real(kind=dp),dimension(:) :: array
    integer(kind=spint) :: i
open(unit=10,file=filename,status='unknown') 
          do i=1,s1
             write(10,*) array(i)
          end do
    close(10)
  end subroutine rwritetofile


end module inputoutput
