module nearfielding
use globalconstants
use neargrouping
implicit none

contains
subroutine createbemdatastruct(t2p,nt,p,np,area,normal,row,col,ncols,potout, &
triafl,nquad,colch,rowch,volch)

implicit none
  ! inputs:
integer(kind=spint),intent(in) :: nt,np,ncols
real(kind=dp),intent(in),dimension(3,np) :: p
integer(kind=spint),intent(in),dimension(3,nt) :: t2p
real(kind=dp),intent(in),dimension(3,nt) :: normal
real(kind=dp),intent(in),dimension(nt) :: area
integer(kind=spint),intent(in),dimension(ncols) :: col,row
integer(kind=spint),intent(in) :: nquad
!internals
integer(kind=spint) :: nquadrhs,itemp,jtemp,ktemp,ktemp2
real(kind=dp),dimension(nquad) :: qwt,qwtrhs
real(kind=dp),dimension(3,nquad) :: qpt,qptrhs

!outputs

integer(kind=spint),intent(out),dimension(nquad*nt) :: rowch,colch
real(kind=dp),intent(out),dimension(3,nquad*nt) :: triafl
real(kind=dp),intent(out),dimension(nquad*nt) :: volch
real(kind=dp),intent(out),dimension(ncols):: potout



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!bem2.Anear=bem.Anear;
!bem2.nquad=bem.nquad
!bem2.epseff=bem.epseff
!bem2.triafl=bem.triafl
!bem2.nhat=bem.nhat
!bem2.vol=bem.vol;
!bem2.Ach=Ach
potout(:)=0.0_dp

call integrationrules2d(qwt,qpt,nquad)
call nearmat2(p,np,t2p,nt,col,row,ncols,normal,area,potout)


  do jtemp=1,nt
  do itemp=1,nquad
ktemp=itemp+(jtemp-1)*nquad
ktemp2=jtemp+(itemp-1)*nt
 triafl(:,ktemp)=p(:,t2p(1,jtemp))*qpt(1,itemp)+ &
                 p(:,t2p(2,jtemp))*qpt(2,itemp)+ & 
                 p(:,t2p(3,jtemp))*qpt(3,itemp)
rowch(ktemp)=ktemp2
colch(ktemp)=jtemp
volch(ktemp)=qwt(itemp)*area(jtemp)

  end do
  end do
call precorr(triafl,nt,volch,nquad,col,row,ncols,normal,potout)
  do jtemp=1,nt
  do itemp=1,nquad
ktemp=jtemp+(itemp-1)*nt
 triafl(:,ktemp)=p(:,t2p(1,jtemp))*qpt(1,itemp)+ &
                 p(:,t2p(2,jtemp))*qpt(2,itemp)+ & 
                 p(:,t2p(3,jtemp))*qpt(3,itemp)
rowch(ktemp)=ktemp
colch(ktemp)=jtemp
volch(ktemp)=qwt(itemp)*area(jtemp)

  end do
  end do

end

subroutine precorr(triafl,nt,qwt,nquad,col,row,ncols,normal,potout)
integer(kind=spint),intent(in) :: nt,nquad,ncols
integer(kind=spint),intent(in),dimension(ncols) :: col,row
real(kind=dp),intent(in),dimension(3,nquad*nt) :: triafl
real(kind=dp),intent(in),dimension(3,nt):: normal
real(kind=dp),intent(in),dimension(nquad*nt) :: qwt
real(kind=dp),intent(inout),dimension(ncols):: potout

!internal variables
real(kind=dp),dimension(3,nquad) :: rs,ro
integer(kind=spint) :: coilid,jquad,iquad,xpre,indtemp,jndtemp
real(kind=dp):: corrpot,radtmp
real(kind=dp),dimension(3):: corrfld,corrtmp
xpre=0

do coilid=1,ncols
if (xpre/=col(coilid)) then
xpre=col(coilid)
!get source points
do iquad=1,nquad
rs(:,iquad)=triafl(:,iquad+(col(coilid)-1)*nquad)
end do
end if
!get observation points
do iquad=1,nquad
ro(:,iquad)=triafl(:,iquad+(row(coilid)-1)*nquad)
end do
!compute correction
!    corrpot=0.0_dp
    corrfld(1)=0.0_dp
    corrfld(2)=0.0_dp
    corrfld(3)=0.0_dp
	do jquad=1,nquad
	do iquad=1,nquad
corrtmp(:)=ro(:,jquad)-rs(:,iquad)
radtmp=sqrt(corrtmp(1)**2+corrtmp(2)**2+corrtmp(3)**2)
if (radtmp>=1d-14) then

indtemp=iquad+nquad*(col(coilid)-1)
jndtemp=jquad+nquad*(row(coilid)-1)
!corrpot=corrpot-qwt(jndtemp)*qwt(indtemp)/radtmp 
radtmp=1.0_dp/radtmp**3
corrfld(1)=corrfld(1)-qwt(jndtemp)*qwt(indtemp)*corrtmp(1)*radtmp
corrfld(2)=corrfld(2)-qwt(jndtemp)*qwt(indtemp)*corrtmp(2)*radtmp
corrfld(3)=corrfld(3)-qwt(jndtemp)*qwt(indtemp)*corrtmp(3)*radtmp
end if
end do
end do
!save correction
potout(coilid)=potout(coilid)+ &
corrfld(1)*normal(1,row(coilid))+ &
corrfld(2)*normal(2,row(coilid))+ &
corrfld(3)*normal(3,row(coilid))

end do
end subroutine

subroutine nearmat2(p,np,t2p,nt,col,row,ncols,normal,area,potout)
implicit none
integer(kind=spint),intent(in) :: nt,np,ncols
integer(kind=spint),intent(in),dimension(ncols) :: col,row
real(kind=dp),intent(in),dimension(3,nt) :: normal
integer(kind=spint),intent(in),dimension(3,nt) :: t2p

real(kind=dp),intent(in),dimension(nt) :: area
real(kind=dp),intent(in),dimension(3,np) :: p
real(kind=dp),intent(inout),dimension(ncols):: potout
!tempdatastruct
real(kind=dp),dimension(16) :: qwt
real(kind=dp),dimension(3,16) :: qpt
real(kind=dp),dimension(:,:),allocatable :: ro

!internal variables
integer(kind=spint) :: coilid,jquad,iquad,iloo,xpre,indtemp,jndtemp
real(kind=dp):: corrpot,radtmp,u3,v3,A2,u0,v0,w0,w02
real(kind=dp),dimension(3):: corrfld,corrtmp,cen,l,nhat,uhat,vhat,s_m,s_p,t0,Rp,Rm,I,boldI,boldk3,beta,K1,R0,ltmp
real(kind=dp),dimension(3,3):: dir,mmmm,q


call integrationrules2d(qwt,qpt,nquadnear)
xpre=0

allocate(ro(3,nquadnear))
do coilid=1,ncols

if (xpre/=col(coilid)) then
xpre=col(coilid)

q(1,1)=p(1,t2p(1,col(coilid)))
q(1,2)=p(2,t2p(1,col(coilid)))
q(1,3)=p(3,t2p(1,col(coilid)))
q(2,1)=p(1,t2p(2,col(coilid)))
q(2,2)=p(2,t2p(2,col(coilid)))
q(2,3)=p(3,t2p(2,col(coilid)))
q(3,1)=p(1,t2p(3,col(coilid)))
q(3,2)=p(2,t2p(3,col(coilid)))
q(3,3)=p(3,t2p(3,col(coilid)))

!geometric processing lines

cen(1)=(q(1,1)+q(2,1)+q(3,1))*0.33333333333333333_dp !centroid
cen(2)=(q(1,2)+q(2,2)+q(3,2))*0.33333333333333333_dp !centroid
cen(3)=(q(1,3)+q(2,3)+q(3,3))*0.33333333333333333_dp !centroid
dir(3,:)=q(2,:)-q(1,:)!construct integration paths
dir(2,:)=q(1,:)-q(3,:)
dir(1,:)=q(3,:)-q(2,:)
l(3)=sqrt(dir(3,1)*dir(3,1)+dir(3,2)*dir(3,2)+dir(3,3)*dir(3,3))!length of paths
l(2)=sqrt(dir(2,1)*dir(2,1)+dir(2,2)*dir(2,2)+dir(2,3)*dir(2,3))!length of paths
l(1)=sqrt(dir(1,1)*dir(1,1)+dir(1,2)*dir(1,2)+dir(1,3)*dir(1,3))!length of paths

! local coordinate system calc (uhat,vhat,nhat)
nhat=(/dir(2,2)*dir(3,3)-dir(2,3)*dir(3,2), &
dir(2,3)*dir(3,1)-dir(2,1)*dir(3,3), &
dir(2,1)*dir(3,2)-dir(2,2)*dir(3,1)/)
A2=sqrt(nhat(1)*nhat(1)+nhat(2)*nhat(2)+nhat(3)*nhat(3))!twicearea
nhat=nhat/A2
ltmp(3)=1/l(3)
uhat(1)=dir(3,1)*ltmp(3)
uhat(2)=dir(3,2)*ltmp(3)
uhat(3)=dir(3,3)*ltmp(3)
vhat=(/nhat(2)*uhat(3)-nhat(3)*uhat(2), &
nhat(3)*uhat(1)-nhat(1)*uhat(3), &
nhat(1)*uhat(2)-nhat(2)*uhat(1)/)
u3=-(dir(2,1)*dir(3,1)+dir(2,2)*dir(3,2)+dir(2,3)*dir(3,3))*ltmp(3)
v3=A2*ltmp(3)

mmmm(:,1)=(q((/2,3,1/),1)-cen(1))
mmmm(:,2)=(q((/2,3,1/),2)-cen(2))
mmmm(:,3)=(q((/2,3,1/),3)-cen(3))
ltmp=(mmmm(:,1)*dir(:,1)+ &
mmmm(:,2)*dir(:,2)+mmmm(:,3)*dir(:,3))/(l*l)
mmmm(1,:)=mmmm(1,:)-ltmp(1)*dir(1,:)
mmmm(1,:)=mmmm(1,:)/sqrt(mmmm(1,1)*mmmm(1,1)+ &
mmmm(1,2)*mmmm(1,2)+mmmm(1,3)*mmmm(1,3))
mmmm(2,:)=mmmm(2,:)-ltmp(2)*dir(2,:)
mmmm(2,:)=mmmm(2,:)/sqrt(mmmm(2,1)*mmmm(2,1)+ &
mmmm(2,2)*mmmm(2,2)+mmmm(2,3)*mmmm(2,3))
mmmm(3,:)=mmmm(3,:)-ltmp(3)*dir(3,:)
mmmm(3,:)=mmmm(3,:)/sqrt(mmmm(3,1)*mmmm(3,1)+ &
mmmm(3,2)*mmmm(3,2)+mmmm(3,3)*mmmm(3,3))


end if

do iquad=1,nquadnear
ro(:,iquad)=p(:,t2p(1,row(coilid)))*qpt(1,iquad)+ &
            p(:,t2p(2,row(coilid)))*qpt(2,iquad)+ & 
            p(:,t2p(3,row(coilid)))*qpt(3,iquad)
end do

 !   corrpot=0.0_dp
    corrfld(1)=0.0_dp
    corrfld(2)=0.0_dp
    corrfld(3)=0.0_dp

boldk3=0.0_dp
ltmp=1.0/l

	do jquad=1,nquadnear

w0=(ro(1,jquad)-q(1,1))*nhat(1)+(ro(2,jquad)-q(1,2))*nhat(2)+(ro(3,jquad)-q(1,3))*nhat(3)
u0=(ro(1,jquad)-q(1,1))*uhat(1)+(ro(2,jquad)-q(1,2))*uhat(2)+(ro(3,jquad)-q(1,3))*uhat(3)
v0=(ro(1,jquad)-q(1,1))*vhat(1)+(ro(2,jquad)-q(1,2))*vhat(2)+(ro(3,jquad)-q(1,3))*vhat(3)
s_m(1)=-((l(3)-u0)*(l(3)-u3)+v0*v3)*ltmp(1)
s_m(2)=-(u3*(u3-u0)+v3*(v3-v0))*ltmp(2)
s_m(3)=-u0
s_p=s_m+l
t0(1)=(v0*(u3-l(3))+v3*(l(3)-u0))*ltmp(1)
t0(2)=(u0*v3-v0*u3)*ltmp(2)
t0(3)=v0
w02=w0*w0
R0(1)=sqrt(t0(1)*t0(1)+w02)
R0(2)=sqrt(t0(2)*t0(2)+w02)
R0(3)=sqrt(t0(3)*t0(3)+w02)
Rp(1)=sqrt((ro(1,jquad)-q(3,1))**2+(ro(2,jquad)-q(3,2))**2+(ro(3,jquad)-q(3,3))**2)
Rp(2)=sqrt((ro(1,jquad)-q(1,1))**2+(ro(2,jquad)-q(1,2))**2+(ro(3,jquad)-q(1,3))**2)
Rp(3)=sqrt((ro(1,jquad)-q(2,1))**2+(ro(2,jquad)-q(2,2))**2+(ro(3,jquad)-q(2,3))**2)
Rm(1)=Rp(3)
Rm(2)=Rp(1)
Rm(3)=Rp(2)
I(1)=(log((Rp(1)+s_p(1))/(Rm(1)+s_m(1))))
I(2)=(log((Rp(2)+s_p(2))/(Rm(2)+s_m(2))))
I(3)=(log((Rp(3)+s_p(3))/(Rm(3)+s_m(3))))
boldI(1)=mmmm(1,1)*I(1)+mmmm(2,1)*I(2)+mmmm(3,1)*I(3)
boldI(2)=mmmm(1,2)*I(1)+mmmm(2,2)*I(2)+mmmm(3,2)*I(3)
boldI(3)=mmmm(1,3)*I(1)+mmmm(2,3)*I(2)+mmmm(3,3)*I(3)

if (w0==0.0_dp) then
    K1(1)=0.0_dp
  !  K1(2)=(t0(1)*I(1)+t0(2)*I(2)+t0(3)*I(3))
else
beta(1)=atan((t0(1)*s_p(1))/(R0(1)**2+abs(w0)*Rp(1)))-atan((t0(1)*s_m(1))/(R0(1)**2+abs(w0)*Rm(1)))
beta(2)=atan((t0(2)*s_p(2))/(R0(2)**2+abs(w0)*Rp(2)))-atan((t0(2)*s_m(2))/(R0(2)**2+abs(w0)*Rm(2)))
beta(3)=atan((t0(3)*s_p(3))/(R0(3)**2+abs(w0)*Rp(3)))-atan((t0(3)*s_m(3))/(R0(3)**2+abs(w0)*Rm(3)))
             K1(1)=1/abs(w0)*(beta(1)+beta(2)+beta(3))
!K1(2)=-(w02)*K1(1)+(t0(1)*I(1)+t0(2)*I(2)+t0(3)*I(3))
end if
boldk3=boldk3-qwt(jquad)*(w0*nhat*K1(1)+boldI)
end do

if (row(coilid)==col(coilid)) then
else
potout(coilid)=potout(coilid)+ &
area(row(coilid))*((-boldk3(1))*normal(1,row(coilid))+ &
                   (-boldk3(2))*normal(2,row(coilid))+ &
                   (-boldk3(3))*normal(3,row(coilid)))
end if
end do
end subroutine



end module nearfielding

