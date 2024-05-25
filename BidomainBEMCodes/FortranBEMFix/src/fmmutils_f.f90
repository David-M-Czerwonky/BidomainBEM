module fmmutils
use globalconstants
implicit none
contains

subroutine computerhs(coil,triamesh,rhs)

type(coilarr),intent(in) :: coil
type(trianglemesh),intent(in) :: triamesh
real(kind=dp),dimension(:),intent(out) :: rhs
real(kind=dp),dimension(:,:),allocatable :: Eprimary,ro

!!!!!!!!!!!internals
real(kind=dp),dimension(3) :: qwt2
real(kind=dp),dimension(3,3) :: qpt2
integer(kind=spint) :: nquad2,i,j,ie



nquad2=3
call integrationrules2d(qwt2,qpt2,nquad2)

allocate(Eprimary(3,triamesh%nt*nquad2))
allocate(ro(3,triamesh%nt*nquad2))
ie=0
do i=1,nquad2
do j=1,triamesh%nt
ro(1:3,ie+j)=triamesh%p(1:3,triamesh%t2p(1,j))*qpt2(1,i)+ &
             triamesh%p(1:3,triamesh%t2p(2,j))*qpt2(2,i)+ &
             triamesh%p(1:3,triamesh%t2p(3,j))*qpt2(3,i)
end do
ie=ie+triamesh%nt
end do
Eprimary=0.0_dp
call computeEprimary(coil%rs,coil%js,ro,Eprimary)
ie=0
rhs=0.0_dp
do i=1,nquad2
do j=1,triamesh%nt
rhs(j)=rhs(j)+qwt2(i)*( &
triamesh%nhat(1,j)*Eprimary(1,ie+j) + &
triamesh%nhat(2,j)*Eprimary(2,ie+j) + &
triamesh%nhat(3,j)*Eprimary(3,ie+j))
end do
ie=ie+triamesh%nt
end do
rhs=rhs*4.0_dp*pi*eps0*triamesh%reg
deallocate(Eprimary)
deallocate(ro)

end subroutine computerhs


subroutine computerhsks(coil,triamesh,rhs)

type(coilarr),intent(in) :: coil
type(trianglemesh),intent(in) :: triamesh
real(kind=dp),dimension(:),intent(out) :: rhs
real(kind=dp),dimension(:,:),allocatable :: Eprimary,ro

!!!!!!!!!!!internals
real(kind=dp),dimension(3) :: qwt2
real(kind=dp),dimension(3,3) :: qpt2
integer(kind=spint) :: nquad2,i,j,ie



nquad2=3
call integrationrules2d(qwt2,qpt2,nquad2)

allocate(Eprimary(3,triamesh%nt*nquad2))
allocate(ro(3,triamesh%nt*nquad2))
ie=0
do i=1,nquad2
do j=1,triamesh%nt
ro(1:3,ie+j)=triamesh%p(1:3,triamesh%t2p(1,j))*qpt2(1,i)+ &
             triamesh%p(1:3,triamesh%t2p(2,j))*qpt2(2,i)+ &
             triamesh%p(1:3,triamesh%t2p(3,j))*qpt2(3,i)
end do
ie=ie+triamesh%nt
end do
Eprimary=0.0_dp
call computeHprimary(coil%rs,coil%js,ro,Eprimary,triamesh%nt,coil%npts,iprecEp)
ie=0
rhs=0.0_dp
do i=1,nquad2
do j=1,triamesh%nt
rhs(j)=rhs(j)+qwt2(i)*( &
triamesh%nhat(1,j)*Eprimary(1,ie+j) + &
triamesh%nhat(2,j)*Eprimary(2,ie+j) + &
triamesh%nhat(3,j)*Eprimary(3,ie+j))
end do
ie=ie+triamesh%nt
end do
rhs=-rhs*4.0_dp*pi*eps0*triamesh%reg
deallocate(Eprimary)
deallocate(ro)

end subroutine computerhsks

subroutine fmmcharge(chspace,triafl,nfmm,fld)
implicit none
real(kind=dp),dimension(:,:),intent(in) :: triafl
real(kind=dp),dimension(:),intent(in) :: chspace
real(kind=dp),dimension(:,:),intent(out) :: fld
integer(kind=spint),intent(in) :: nfmm
real(kind=dp),dimension(nfmm) :: pot

          call lfmm3d_s_c_g(iprecEP,nfmm,triafl,chspace,pot,fld)

fld=-fld
end subroutine fmmcharge


subroutine fmmchargetarget(chspace,triafl,nfmm,rtarg,fldtarg,ntarget)
implicit none
integer(kind=spint),intent(in) :: nfmm,ntarget
real(kind=dp),dimension(3,nfmm),intent(in) :: triafl
real(kind=dp),dimension(nfmm),intent(in) :: chspace
real(kind=dp),dimension(3,ntarget),intent(in) :: rtarg
real(kind=dp),dimension(3,ntarget),intent(inout) :: fldtarg
real(kind=dp),dimension(ntarget) :: pottarg

call lfmm3d_t_c_g(iprecEP,nfmm,triafl,chspace,ntarget,rtarg,pottarg,fldtarg)
fldtarg=-fldtarg
end subroutine fmmchargetarget


subroutine computeEprimary(rs,js,robs,Eprimary)
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: spint=kind(8)
  real(kind=dp),parameter :: pi=3.141592653589793_dp,epsthresh=1d-12
  real(kind=dp),parameter ::eps0=8.854187820000000d-12
  real(kind=dp) ::  mu0over4pi=-1.000000000000000d-7
real(kind=dp),dimension(:,:),intent(in) :: rs
real(kind=dp),dimension(:,:),intent(in) :: robs
real(kind=dp),dimension(:,:),intent(in) :: js
real(kind=dp),dimension(:,:),intent(out) :: Eprimary
integer :: nsource,ntarget

integer(kind=spint) :: i,j,ier

nsource=size(rs,2)
ntarget=size(robs,2)

call lfmm3d_t_c_p_vec(3,iprecEP,nsource,rs,js,ntarget,robs,Eprimary)
Eprimary=mu0over4pi*Eprimary



end subroutine


subroutine computeHprimary(rs,js,robs,Hprimary,ntarget,nsource,precEP)
  use omp_lib
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: spint=kind(8)
  real(kind=dp) ::  mu0over4pi
real(kind=dp),dimension(3,nsource),intent(in) :: rs
real(kind=dp),dimension(3,ntarget),intent(in) :: robs
real(kind=dp),dimension(3,nsource),intent(in) :: js
real(kind=dp),dimension(3,ntarget),intent(out) :: Hprimary
integer :: nsource,ntarget
real(kind=dp),intent(in) :: precEP

integer(kind=spint) :: i,j,ier
real(kind=dp),dimension(:),allocatable :: pottarg,tempcharge
real(kind=dp),dimension(:,:),allocatable :: fldtarg1

  !f2py threadsafe
  !f2py intent(hide) :: nsource=shape(rs, 1), ntarget=shape(robs, 1)
  !f2py depend(ntarget) Hprimary
mu0over4pi=-1.000000000000000d-7

allocate(tempcharge(nsource))
allocate(fldtarg1(3,ntarget))
allocate(pottarg(ntarget))


tempcharge=reshape(mu0over4pi*js(1,:),(/nsource/))
call lfmm3d_t_c_g(precEP,nsource,rs,tempcharge,ntarget,robs,pottarg,fldtarg1)


	Hprimary(2,:)=-reshape(fldtarg1(3,:),(/ntarget/))
	Hprimary(3,:)=reshape(fldtarg1(2,:),(/ntarget/))

tempcharge=reshape(mu0over4pi*js(2,:),(/nsource/))
call lfmm3d_t_c_g(precEP,nsource,rs,tempcharge,ntarget,robs,pottarg,fldtarg1)

	Hprimary(1,:)=reshape(fldtarg1(3,:),(/ntarget/))
	Hprimary(3,:)=Hprimary(3,:)-reshape(fldtarg1(1,:),(/ntarget/))

tempcharge=reshape(mu0over4pi*js(3,:),(/nsource/))
call lfmm3d_t_c_g(precEP,nsource,rs,tempcharge,ntarget,robs,pottarg,fldtarg1)

	Hprimary(1,:)=Hprimary(1,:)-reshape(fldtarg1(2,:),(/ntarget/))
	Hprimary(2,:)=Hprimary(2,:)+reshape(fldtarg1(1,:),(/ntarget/))


deallocate(fldtarg1)
deallocate(pottarg)
deallocate(tempcharge)
end subroutine



end module fmmutils
