module globalconstants
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: spint=kind(8)
  real(kind=dp),parameter :: pi=3.141592653589793_dp,epsthresh=1d-12
  real(kind=dp),parameter ::eps0=8.854187820000000d-12
  real(kind=dp) ::  mu0=1.256637060000000d-6
  real(kind=dp) :: iprecEp=1.0d-3,iprecmatvec=1.0d-3
  integer(kind=spint) :: nquad=3,nquadnear=16
type coilarr
!curvec has to be rescaled with mu0/(4*pi)
integer(kind=spint) :: npts
real(kind=dp), dimension(:,:),allocatable :: rs
real(kind=dp), dimension(:,:),allocatable :: js
end type coilarr

type trianglemesh
integer(kind=spint) :: np,nt
real(kind=dp), dimension(:,:),allocatable :: p,nhat
real(kind=dp), dimension(:),allocatable :: reg,vol
integer(kind=spint), dimension(:,:),allocatable :: t2p
end type trianglemesh

type sparse
integer(kind=spint) :: npts
real(kind=dp), dimension(:),allocatable :: val
integer(kind=spint), dimension(:),allocatable :: col,row
integer(kind=spint) :: nnz,nrows,ncols
end type sparse

contains

subroutine integrationrules2d(qwt,qpt,nquadi)
implicit none
integer(kind=spint),intent(in) :: nquadi
real(kind=dp),dimension(nquadi),intent(out) :: qwt
real(kind=dp),dimension(3,nquadi),intent(out) :: qpt


if (nquadi==1) then
qwt=1.0_dp
qpt(:,1)=(/ 1.0_dp/3.0_dp,1.0_dp/3.0_dp,1.0_dp/3.0_dp /)
elseif (nquadi==3)then

qwt=(/ 1.0_dp/3.0_dp,1.0_dp/3.0_dp,1.0_dp/3.0_dp /)
qpt(:,1)=(/ 1.0_dp/6.0_dp,1.0_dp/6.0_dp,2.0_dp/3.0_dp /)
qpt(:,2)=(/ 1.0_dp/6.0_dp,2.0_dp/3.0_dp,1.0_dp/6.0_dp /)
qpt(:,3)=(/ 2.0_dp/3.0_dp,1.0_dp/6.0_dp,1.0_dp/6.0_dp /)

elseif (nquadi==16) then

qwt=(/ &
0.1443156076777900_dp,0.0950916342672800_dp,0.0950916342672800_dp,0.09509163426728000_dp, &
0.1032173705347200_dp,0.1032173705347200_dp,0.1032173705347200_dp,0.03245849762320000_dp, &
0.0324584976232000_dp,0.0324584976232000_dp,0.0272303141744300_dp,0.02723031417443000_dp, &
0.0272303141744300_dp,0.0272303141744300_dp,0.0272303141744300_dp,0.02723031417443000_dp /)
qpt(:,1)=(/0.333333333333330_dp,0.333333333333330_dp,0.333333333333340_dp /)
qpt(:,2)=(/0.459292588292720_dp,0.459292588292720_dp,0.081414823414560_dp /)
qpt(:,3)=(/0.459292588292720_dp,0.081414823414550_dp,0.459292588292730_dp /)
qpt(:,4)=(/0.081414823414550_dp,0.459292588292720_dp,0.459292588292730_dp /)
qpt(:,5)=(/0.170569307751760_dp,0.170569307751760_dp,0.658861384496480_dp /)
qpt(:,6)=(/0.170569307751760_dp,0.658861384496480_dp,0.170569307751760_dp /)
qpt(:,7)=(/0.658861384496480_dp,0.170569307751760_dp,0.170569307751760_dp /)
qpt(:,8)=(/0.050547228317030_dp,0.050547228317030_dp,0.898905543365940_dp /)
qpt(:,9)=(/0.050547228317030_dp,0.898905543365940_dp,0.050547228317030_dp /)
qpt(:,10)=(/0.898905543365940_dp,0.050547228317030_dp,0.050547228317030_dp /)
qpt(:,11)=(/0.263112829634640_dp,0.728492392955400_dp,0.008394777409960_dp /)
qpt(:,12)=(/0.728492392955400_dp,0.008394777409960_dp,0.263112829634640_dp /)
qpt(:,13)=(/0.008394777409960_dp,0.263112829634640_dp,0.728492392955400_dp /)
qpt(:,14)=(/0.728492392955400_dp,0.263112829634640_dp,0.008394777409960_dp /)
qpt(:,15)=(/0.263112829634640_dp,0.008394777409960_dp,0.728492392955400_dp /)
qpt(:,16)=(/0.008394777409960_dp,0.728492392955400_dp,0.263112829634640_dp /)
end if
end subroutine integrationrules2d

  function isnan (x)
    logical :: isnan
    real(kind=dp) ::  x
    if (x.ne.x) then
       isnan=.TRUE.
    else
       isnan=.FALSE.
    end if
    return
  end function isnan
end module globalconstants
