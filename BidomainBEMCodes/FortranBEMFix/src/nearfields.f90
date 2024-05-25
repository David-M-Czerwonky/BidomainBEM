module nearfielding
use globalconstants
use neargrouping
!use fmmutils
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
real(kind=dp),intent(out),dimension(ncols) :: potout



volch(:)=0.0_dp
potout(:)=0.0_dp


call integrationrules2d(qwt,qpt,nquad)
!call nearmat2(p,np,t2p,nt,col,row,ncols,normal,area,potout)
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

subroutine nearmat2(p,np,t2p,nt,col,row,ncols,normal,area,Dout,Sout,Doutlin,Soutlin)
implicit none
integer(kind=spint),intent(in) :: nt,np,ncols
integer(kind=spint),intent(in),dimension(ncols) :: col,row
real(kind=dp),intent(in),dimension(3,nt) :: normal
integer(kind=spint),intent(in),dimension(3,nt) :: t2p

real(kind=dp),intent(in),dimension(nt) :: area
real(kind=dp),intent(in),dimension(3,np) :: p
real(kind=dp),intent(out),dimension(ncols):: Sout,Dout
real(kind=dp),intent(out),dimension(ncols,9):: Soutlin, Doutlin
!tempdatastruct
real(kind=dp),dimension(16) :: qwt
real(kind=dp),dimension(3,16) :: qpt
real(kind=dp),dimension(:,:),allocatable :: ro

!internal variables
integer(kind=spint) :: coilid,jquad,iquad,iloo,xpre,indtemp,jndtemp
real(kind=dp):: corrpot,radtmp,u3,v3,A2,u0,v0,w0,I1,Iu,Iv,N10,N20,N30,f1t,f2t,f3t
real(kind=dp),dimension(3):: corrfld,corrtmp,cen,l,nhat,uhat,vhat,s_m,s_p,t0,Rp,Rm,I,boldI,boldk3,beta,K1,R0,f3i,bk3,boldUa,boldVa
real(kind=dp),dimension(3,3):: dir,mmmm,q,f_i,f1,f2



call integrationrules2d(qwt,qpt,nquadnear)
xpre=0

allocate(ro(3,nquadnear))
do coilid=1,ncols

if (xpre/=col(coilid)) then
xpre=col(coilid)
do iquad=1,3
q(iquad,1)=p(1,t2p(iquad,col(coilid)))
q(iquad,2)=p(2,t2p(iquad,col(coilid)))
q(iquad,3)=p(3,t2p(iquad,col(coilid)))
end do

!geometric processing lines

cen(1)=(q(1,1)+q(2,1)+q(3,1))/3.0_dp !centroid
cen(2)=(q(1,2)+q(2,2)+q(3,2))/3.0_dp !centroid
cen(3)=(q(1,3)+q(2,3)+q(3,3))/3.0_dp !centroid
dir(3,:)=q(2,:)-q(1,:)!construct integration paths
dir(2,:)=q(1,:)-q(3,:)
dir(1,:)=q(3,:)-q(2,:)
l(3)=sqrt(dir(3,1)**2.0_dp+dir(3,2)**2.0_dp+dir(3,3)**2.0_dp)!length of paths
l(2)=sqrt(dir(2,1)**2.0_dp+dir(2,2)**2.0_dp+dir(2,3)**2.0_dp)!length of paths
l(1)=sqrt(dir(1,1)**2.0_dp+dir(1,2)**2.0_dp+dir(1,3)**2.0_dp)!length of paths
! local coordinate system calc (uhat,vhat,nhat)
nhat=(/dir(2,2)*dir(3,3)-dir(2,3)*dir(3,2), &
dir(2,3)*dir(3,1)-dir(2,1)*dir(3,3), &
dir(2,1)*dir(3,2)-dir(2,2)*dir(3,1)/)
A2=sqrt(nhat(1)**2.0_dp+nhat(2)**2.0_dp+nhat(3)**2.0_dp)!twicearea
nhat=nhat/A2
uhat(1)=dir(3,1)/l(3)
uhat(2)=dir(3,2)/l(3)
uhat(3)=dir(3,3)/l(3)
vhat=(/nhat(2)*uhat(3)-nhat(3)*uhat(2), &
nhat(3)*uhat(1)-nhat(1)*uhat(3), &
nhat(1)*uhat(2)-nhat(2)*uhat(1)/)
u3=-(dir(2,1)*dir(3,1)+dir(2,2)*dir(3,2)+dir(2,3)*dir(3,3))/l(3)
v3=A2/l(3)
do iloo=1,3
mmmm(iloo,1)=(q(mod(iloo,3)+1,1)-cen(1))
mmmm(iloo,2)=(q(mod(iloo,3)+1,2)-cen(2))
mmmm(iloo,3)=(q(mod(iloo,3)+1,3)-cen(3))
dir(iloo,:)=dir(iloo,:)/l(iloo)
mmmm(iloo,:)=mmmm(iloo,:)-(mmmm(iloo,1)*dir(iloo,1)+ &
mmmm(iloo,2)*dir(iloo,2)+mmmm(iloo,3)*dir(iloo,3))*dir(iloo,:)
mmmm(iloo,:)=mmmm(iloo,:)/sqrt(mmmm(iloo,1)**2.0_dp+ &
mmmm(iloo,2)**2.0_dp+mmmm(iloo,3)**2.0_dp);
end do
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
I1=0.0_dp
boldk3=0.0_dp
f1=0.0_dp
f2=0.0_dp

	do jquad=1,nquadnear
u0=(ro(1,jquad)-q(1,1))*uhat(1)+(ro(2,jquad)-q(1,2))*uhat(2)+(ro(3,jquad)-q(1,3))*uhat(3)
v0=(ro(1,jquad)-q(1,1))*vhat(1)+(ro(2,jquad)-q(1,2))*vhat(2)+(ro(3,jquad)-q(1,3))*vhat(3)
w0=(ro(1,jquad)-q(1,1))*nhat(1)+(ro(2,jquad)-q(1,2))*nhat(2)+(ro(3,jquad)-q(1,3))*nhat(3)
s_m(1)=-((l(3)-u0)*(l(3)-u3)+v0*v3)/l(1)
s_m(2)=-(u3*(u3-u0)+v3*(v3-v0))/l(2)
s_m(3)=-u0
s_p=s_m+l
t0(1)=(v0*(u3-l(3))+v3*(l(3)-u0))/l(1)
t0(2)=(u0*v3-v0*u3)/l(2)
t0(3)=v0
R0(1)=sqrt(t0(1)**2.0_dp+w0**2.0_dp)
R0(2)=sqrt(t0(2)**2.0_dp+w0**2.0_dp)
R0(3)=sqrt(t0(3)**2.0_dp+w0**2.0_dp)
Rp(1)=sqrt((ro(1,jquad)-q(3,1))**2.0_dp+(ro(2,jquad)-q(3,2))**2.0_dp+(ro(3,jquad)-q(3,3))**2.0_dp)
Rp(2)=sqrt((ro(1,jquad)-q(1,1))**2.0_dp+(ro(2,jquad)-q(1,2))**2.0_dp+(ro(3,jquad)-q(1,3))**2.0_dp)
Rp(3)=sqrt((ro(1,jquad)-q(2,1))**2.0_dp+(ro(2,jquad)-q(2,2))**2.0_dp+(ro(3,jquad)-q(2,3))**2.0_dp)
Rm(1)=Rp(3)
Rm(2)=Rp(1)
Rm(3)=Rp(2)
I(1)=(log((Rp(1)+s_p(1))/(Rm(1)+s_m(1))))
I(2)=(log((Rp(2)+s_p(2))/(Rm(2)+s_m(2))))
I(3)=(log((Rp(3)+s_p(3))/(Rm(3)+s_m(3))))
if (I(1)/=I(1)) then
I(1)=0.0_dp;
else if (I(1)-1.0_dp==I(1)) then
I(1)=0.0_dp
end if
if (I(2)/=I(2)) then
I(2)=0.0_dp;
else if (I(2)-1.0_dp==I(2)) then
I(2)=0.0_dp
end if
if (I(3)/=I(3)) then
I(3)=0.0_dp
else if (I(3)-1.0_dp==I(3)) then
I(3)=0.0_dp
end if
boldI(1)=mmmm(1,1)*I(1)+mmmm(2,1)*I(2)+mmmm(3,1)*I(3)
boldI(2)=mmmm(1,2)*I(1)+mmmm(2,2)*I(2)+mmmm(3,2)*I(3)
boldI(3)=mmmm(1,3)*I(1)+mmmm(2,3)*I(2)+mmmm(3,3)*I(3)

f3i(3)=s_p(3)*Rp(3)-s_m(3)*Rm(3)+R0(3)**2*I(3)
f3i(1)=s_p(1)*Rp(1)-s_m(1)*Rm(1)+R0(1)**2*I(1)
f3i(2)=s_p(2)*Rp(2)-s_m(2)*Rm(2)+R0(2)**2*I(2)
if (abs(w0)<=1d-13) then
    K1(1)=0.0_dp
    K1(2)=(t0(1)*I(1)+t0(2)*I(2)+t0(3)*I(3))
    do iloo=1,3
    I(iloo)=0.0_dp;
    if ((abs(Rp(iloo)+s_p(iloo))>=1d-12).AND.(abs(Rm(iloo)+s_m(iloo))>=1d-12)) then
    I(iloo)=(log((Rp(iloo)+s_p(iloo))/(Rm(iloo)+s_m(iloo))))
    end if
	if (abs(t0(iloo))<=1d-13) then
	beta(iloo)=0.0_dp
	end if
    end do
    boldI(1)=mmmm(1,1)*I(1)+mmmm(2,1)*I(2)+mmmm(3,1)*I(3)
    boldI(2)=mmmm(1,2)*I(1)+mmmm(2,2)*I(2)+mmmm(3,2)*I(3)
    boldI(3)=mmmm(1,3)*I(1)+mmmm(2,3)*I(2)+mmmm(3,3)*I(3)


else
beta(1)=atan((t0(1)*s_p(1))/(R0(1)**2.0_dp+abs(w0)*Rp(1)))-atan((t0(1)*s_m(1))/(R0(1)**2.0_dp+abs(w0)*Rm(1)))
beta(2)=atan((t0(2)*s_p(2))/(R0(2)**2.0_dp+abs(w0)*Rp(2)))-atan((t0(2)*s_m(2))/(R0(2)**2.0_dp+abs(w0)*Rm(2)))
beta(3)=atan((t0(3)*s_p(3))/(R0(3)**2.0_dp+abs(w0)*Rp(3)))-atan((t0(3)*s_m(3))/(R0(3)**2.0_dp+abs(w0)*Rm(3)))
    do iloo=1,3
	if (abs(t0(iloo))<=1d-13) then
	beta(iloo)=0.0_dp
	end if
    end do
             K1(1)=1/abs(w0)*(beta(1)+beta(2)+beta(3))
K1(2)=-(w0**2.0_dp)*K1(1)+(t0(1)*I(1)+t0(2)*I(2)+t0(3)*I(3))
end if

bk3=-(w0*nhat*K1(1)+boldI)
Iu=0.5_dp*((mmmm(1,1)*uhat(1)+mmmm(1,2)*uhat(2)+mmmm(1,3)*uhat(3))*f3i(1)+&
	       (mmmm(2,1)*uhat(1)+mmmm(2,2)*uhat(2)+mmmm(2,3)*uhat(3))*f3i(2)+&
           (mmmm(3,1)*uhat(1)+mmmm(3,2)*uhat(2)+mmmm(3,3)*uhat(3))*f3i(3))

Iv=0.5_dp*((mmmm(1,1)*vhat(1)+mmmm(1,2)*vhat(2)+mmmm(1,3)*vhat(3))*f3i(1)+&
	       (mmmm(2,1)*vhat(1)+mmmm(2,2)*vhat(2)+mmmm(2,3)*vhat(3))*f3i(2)+&
           (mmmm(3,1)*vhat(1)+mmmm(3,2)*vhat(2)+mmmm(3,3)*vhat(3))*f3i(3))
f_i(1,:)=dir(1,:)*t0(1)*I(1)-mmmm(1,:)*(Rp(1)-Rm(1))
f_i(2,:)=dir(2,:)*t0(2)*I(2)-mmmm(2,:)*(Rp(2)-Rm(2))
f_i(3,:)=dir(3,:)*t0(3)*I(3)-mmmm(3,:)*(Rp(3)-Rm(3))


boldUa= &
(dir(1,1)*uhat(1)+dir(1,2)*uhat(2)+dir(1,3)*uhat(3))*f_i(1,:)+ &
(dir(2,1)*uhat(1)+dir(2,2)*uhat(2)+dir(2,3)*uhat(3))*f_i(2,:)+ &
(dir(3,1)*uhat(1)+dir(3,2)*uhat(2)+dir(3,3)*uhat(3))*f_i(3,:)+ &
nhat*w0*(boldI(1)*uhat(1)+boldI(2)*uhat(2)+boldI(3)*uhat(3))- &
uhat*abs(w0)*sum(beta)

boldVa= &
(dir(1,1)*vhat(1)+dir(1,2)*vhat(2)+dir(1,3)*vhat(3))*f_i(1,:)+ &
(dir(2,1)*vhat(1)+dir(2,2)*vhat(2)+dir(2,3)*vhat(3))*f_i(2,:)+ &
(dir(3,1)*vhat(1)+dir(3,2)*vhat(2)+dir(3,3)*vhat(3))*f_i(3,:)+ &
nhat*w0*(boldI(1)*vhat(1)+boldI(2)*vhat(2)+boldI(3)*vhat(3))- &
vhat*abs(w0)*sum(beta)

N10=1-u0/l(3)+(u3/l(3)-1)*v0/v3
N20=u0/l(3)-(u3/l(3))*v0/v3
N30=v0/v3


f1(:,1)=f1(:,1)+qwt(jquad)*qpt(:,jquad)* &
((N10*bk3(1)-boldUa(1)/l(3)+(u3/l(3)-1)*boldVa(1)/v3)*normal(1,row(coilid))+ &
 (N10*bk3(2)-boldUa(2)/l(3)+(u3/l(3)-1)*boldVa(2)/v3)*normal(2,row(coilid))+ &
 (N10*bk3(3)-boldUa(3)/l(3)+(u3/l(3)-1)*boldVa(3)/v3)*normal(3,row(coilid)))
f1(:,2)=f1(:,2)+qwt(jquad)*qpt(:,jquad)* &
((N20*bk3(1)+boldUa(1)/l(3)-u3/l(3)*boldVa(1)/v3)*normal(1,row(coilid))+ &
 (N20*bk3(2)+boldUa(2)/l(3)-u3/l(3)*boldVa(2)/v3)*normal(2,row(coilid))+ &
 (N20*bk3(3)+boldUa(3)/l(3)-u3/l(3)*boldVa(3)/v3)*normal(3,row(coilid)))
f1(:,3)=f1(:,3)+qwt(jquad)*qpt(:,jquad)* &
((N30*bk3(1)+boldVa(1)/v3)*normal(1,row(coilid))+ &
 (N30*bk3(2)+boldVa(2)/v3)*normal(2,row(coilid))+ &
 (N30*bk3(3)+boldVa(3)/v3)*normal(3,row(coilid)))

f2(:,1)=f2(:,1)+qwt(jquad)*qpt(:,jquad)*(N10*K1(2)-Iu/l(3)+(u3/l(3)-1)*Iv/v3)
f2(:,2)=f2(:,2)+qwt(jquad)*qpt(:,jquad)*(N20*K1(2)+Iu/l(3)-u3/l(3)*Iv/v3)
f2(:,3)=f2(:,3)+qwt(jquad)*qpt(:,jquad)*(N30*K1(2)+Iv/v3)


boldk3=boldk3+qwt(jquad)*bk3
I1=I1+qwt(jquad)*K1(2)

end do



if (row(coilid)==col(coilid)) then
  Dout(coilid)=0.0_dp
  Doutlin(coilid,:)=0.0_dp
else
Dout(coilid)=Dout(coilid)+ &
area(row(coilid))*((-boldk3(1))*normal(1,row(coilid))+ &
                   (-boldk3(2))*normal(2,row(coilid))+ &
                   (-boldk3(3))*normal(3,row(coilid)))
Doutlin(coilid,1)=Doutlin(coilid,1)+area(row(coilid))*f1(1,1)
Doutlin(coilid,2)=Doutlin(coilid,2)+area(row(coilid))*f1(2,1)
Doutlin(coilid,3)=Doutlin(coilid,3)+area(row(coilid))*f1(3,1)
Doutlin(coilid,4)=Doutlin(coilid,4)+area(row(coilid))*f1(1,2)
Doutlin(coilid,5)=Doutlin(coilid,5)+area(row(coilid))*f1(2,2)
Doutlin(coilid,6)=Doutlin(coilid,6)+area(row(coilid))*f1(3,2)
Doutlin(coilid,7)=Doutlin(coilid,7)+area(row(coilid))*f1(1,3)
Doutlin(coilid,8)=Doutlin(coilid,8)+area(row(coilid))*f1(2,3)
Doutlin(coilid,9)=Doutlin(coilid,9)+area(row(coilid))*f1(3,3)

end if
Sout(coilid)=Sout(coilid)+area(row(coilid))*I1
Soutlin(coilid,1)=Soutlin(coilid,1)+area(row(coilid))*f2(1,1)
Soutlin(coilid,2)=Soutlin(coilid,2)+area(row(coilid))*f2(2,1)
Soutlin(coilid,3)=Soutlin(coilid,3)+area(row(coilid))*f2(3,1)
Soutlin(coilid,4)=Soutlin(coilid,4)+area(row(coilid))*f2(1,2)
Soutlin(coilid,5)=Soutlin(coilid,5)+area(row(coilid))*f2(2,2)
Soutlin(coilid,6)=Soutlin(coilid,6)+area(row(coilid))*f2(3,2)
Soutlin(coilid,7)=Soutlin(coilid,7)+area(row(coilid))*f2(1,3)
Soutlin(coilid,8)=Soutlin(coilid,8)+area(row(coilid))*f2(2,3)
Soutlin(coilid,9)=Soutlin(coilid,9)+area(row(coilid))*f2(3,3)


end do


end subroutine 



subroutine computeEsecondary(t2p,nt,p,np,area,normal,col,row,ncols,Efldcorr,robs,nobs, &
nquad,xval)

implicit none
  ! inputs:
integer(kind=spint),intent(in) :: nt,np,ncols,nobs
real(kind=dp),intent(in),dimension(3,np) :: p
real(kind=dp),intent(in),dimension(3,nobs) :: robs
integer(kind=spint),intent(in),dimension(3,nt) :: t2p
real(kind=dp),intent(in),dimension(3,nt) :: normal
real(kind=dp),intent(in),dimension(nt) :: area,xval
integer(kind=spint),intent(in),dimension(ncols) :: col,row
integer(kind=spint),intent(in) :: nquad
!internals
integer(kind=spint) :: nquadrhs,itemp,jtemp,ktemp
real(kind=dp),dimension(nquad) :: qwt,qwtrhs
real(kind=dp),dimension(3,nquad) :: qpt,qptrhs
real(kind=dp),dimension(nquad*nt) :: chspace
real(kind=dp),dimension(3,nobs) :: fldvv
!outputs

integer(kind=spint),dimension(nquad*nt) :: rowch,colch
real(kind=dp),dimension(:,:),allocatable :: triafl
real(kind=dp),dimension(nquad*nt) :: volch
real(kind=dp),intent(inout),dimension(3,nobs):: Efldcorr
real(kind=dp),dimension(nt) :: xvaltemp

allocate(triafl(3,nquad*nt))
xvaltemp=xval/(4*pi*eps0)
!!!get integration rules
call integrationrules2d(qwt,qpt,nquad)
!generate quadrature points and weights
  do itemp=1,nquad
  do jtemp=1,nt
ktemp=itemp+(jtemp-1)*nquad
 triafl(:,ktemp)=p(:,t2p(1,jtemp))*qpt(1,itemp)+ &
                 p(:,t2p(2,jtemp))*qpt(2,itemp)+ &
                 p(:,t2p(3,jtemp))*qpt(3,itemp)
rowch(ktemp)=ktemp
colch(ktemp)=jtemp
volch(ktemp)=qwt(itemp)*area(jtemp)
chspace(ktemp)=volch(ktemp)*xvaltemp(jtemp)
  end do
  end do
!generate dipole representation of charges
jtemp=nquad*nt

!use FMM for particle field evaluation

fldvv=0.0_dp
!call fmmchargetarget(chspace,triafl,jtemp,robs,fldvv,nobs)
!Efldcorr=Efldcorr+fldvv
!delete near entries
call precorrfld(triafl,nt,xvaltemp,chspace,nquad,col,row,ncols,normal,Efldcorr,robs,nobs)
!add analytic expression
call nearfld(p,np,t2p,nt,xvaltemp,col,row,ncols,normal,area,Efldcorr,robs,nobs)

end subroutine computeEsecondary

subroutine precorrfld(triafl,nt,xval,qwt,nquad,col,row,ncols,normal,Efldcorr,robs,nobs)
integer(kind=spint),intent(in) :: nt,nquad,ncols,nobs
integer(kind=spint),intent(in),dimension(ncols) :: col,row
real(kind=dp),intent(in),dimension(nt) ::  xval
real(kind=dp),intent(in),dimension(3,nquad*nt) :: triafl
real(kind=dp),intent(in),dimension(3,nt):: normal
real(kind=dp),intent(in),dimension(nquad*nt) :: qwt
real(kind=dp),intent(in),dimension(3,nobs):: robs

real(kind=dp),intent(inout),dimension(3,nobs):: Efldcorr

!internal variables
real(kind=dp),dimension(3,nquad) :: rs
integer(kind=spint) :: coilid,jquad,iquad,xpre,indtemp,jndtemp
real(kind=dp):: corrpot,radtmp
real(kind=dp),dimension(3):: corrfld,corrtmp,ro
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

!compute correction
!    corrpot=0.0_dp
    corrfld(1)=0.0_dp
    corrfld(2)=0.0_dp
    corrfld(3)=0.0_dp
	do iquad=1,nquad
corrtmp(:)=robs(:,row(coilid))-rs(:,iquad)
radtmp=sqrt(corrtmp(1)**2.0_dp+corrtmp(2)**2.0_dp+corrtmp(3)**2.0_dp)
if (radtmp>=1d-14) then

indtemp=iquad+nquad*(col(coilid)-1)
!corrpot=corrpot-qwt(jndtemp)*qwt(indtemp)/radtmp
radtmp=radtmp**(-3.0_dp)
corrfld(1)=corrfld(1)-real(qwt(indtemp))*corrtmp(1)*radtmp
corrfld(2)=corrfld(2)-real(qwt(indtemp))*corrtmp(2)*radtmp
corrfld(3)=corrfld(3)-real(qwt(indtemp))*corrtmp(3)*radtmp
end if
end do
!save correction
Efldcorr(1,row(coilid))=Efldcorr(1,row(coilid))+corrfld(1)
Efldcorr(2,row(coilid))=Efldcorr(2,row(coilid))+corrfld(2)
Efldcorr(3,row(coilid))=Efldcorr(3,row(coilid))+corrfld(3)

end do
end subroutine

subroutine nearfld(p,np,t2p,nt,xval,col,row,ncols,normal,area,Efldcorr,robs,nobs)
implicit none
integer(kind=spint),intent(in) :: nt,np,ncols,nobs
integer(kind=spint),intent(in),dimension(ncols) :: col,row
real(kind=dp),intent(in),dimension(3,nt) :: normal
integer(kind=spint),intent(in),dimension(3,nt) :: t2p

real(kind=dp),intent(in),dimension(nt) :: area,xval
real(kind=dp),intent(in),dimension(3,np) :: p
real(kind=dp),intent(in),dimension(3,nobs):: robs
real(kind=dp),intent(inout),dimension(3,nobs):: Efldcorr
!tempdatastruct
real(kind=dp),dimension(16) :: qwt
real(kind=dp),dimension(3,16) :: qpt
real(kind=dp),dimension(3) :: ro

!internal variables
integer(kind=spint) :: coilid,jquad,iquad,iloo,xpre,indtemp,jndtemp
real(kind=dp):: corrpot,radtmp,u3,v3,A2,u0,v0,w0
real(kind=dp),dimension(3):: corrfld,corrtmp,cen,l,nhat,uhat,vhat,s_m,s_p,t0,Rp,Rm,I,boldI,boldk3,beta,K1,R0
real(kind=dp),dimension(3,3):: dir,mmmm,q


call integrationrules2d(qwt,qpt,nquadnear)
xpre=0

do coilid=1,ncols

if (xpre/=col(coilid)) then
xpre=col(coilid)
do iquad=1,3
q(iquad,1)=p(1,t2p(iquad,col(coilid)))
q(iquad,2)=p(2,t2p(iquad,col(coilid)))
q(iquad,3)=p(3,t2p(iquad,col(coilid)))
end do

!geometric processing lines

cen(1)=(q(1,1)+q(2,1)+q(3,1))/3.0_dp !centroid
cen(2)=(q(1,2)+q(2,2)+q(3,2))/3.0_dp !centroid
cen(3)=(q(1,3)+q(2,3)+q(3,3))/3.0_dp !centroid
dir(3,:)=q(2,:)-q(1,:)!construct integration paths
dir(2,:)=q(1,:)-q(3,:)
dir(1,:)=q(3,:)-q(2,:)
l(3)=sqrt(dir(3,1)**2.0_dp+dir(3,2)**2.0_dp+dir(3,3)**2.0_dp)!length of paths
l(2)=sqrt(dir(2,1)**2.0_dp+dir(2,2)**2.0_dp+dir(2,3)**2.0_dp)!length of paths
l(1)=sqrt(dir(1,1)**2.0_dp+dir(1,2)**2.0_dp+dir(1,3)**2.0_dp)!length of paths
! local coordinate system calc (uhat,vhat,nhat)
nhat=(/dir(2,2)*dir(3,3)-dir(2,3)*dir(3,2), &
dir(2,3)*dir(3,1)-dir(2,1)*dir(3,3), &
dir(2,1)*dir(3,2)-dir(2,2)*dir(3,1)/)
A2=sqrt(nhat(1)**2.0_dp+nhat(2)**2.0_dp+nhat(3)**2.0_dp)!twicearea
nhat=nhat/A2
uhat(1)=dir(3,1)/l(3)
uhat(2)=dir(3,2)/l(3)
uhat(3)=dir(3,3)/l(3)
vhat=(/nhat(2)*uhat(3)-nhat(3)*uhat(2), &
nhat(3)*uhat(1)-nhat(1)*uhat(3), &
nhat(1)*uhat(2)-nhat(2)*uhat(1)/)
u3=-(dir(2,1)*dir(3,1)+dir(2,2)*dir(3,2)+dir(2,3)*dir(3,3))/l(3)
v3=A2/l(3)
do iloo=1,3
mmmm(iloo,1)=(q(mod(iloo,3)+1,1)-cen(1))
mmmm(iloo,2)=(q(mod(iloo,3)+1,2)-cen(2))
mmmm(iloo,3)=(q(mod(iloo,3)+1,3)-cen(3))
mmmm(iloo,:)=mmmm(iloo,:)-(mmmm(iloo,1)*dir(iloo,1)+ &
mmmm(iloo,2)*dir(iloo,2)+mmmm(iloo,3)*dir(iloo,3))*dir(iloo,:)/(l(iloo)**2.0_dp);
mmmm(iloo,:)=mmmm(iloo,:)/sqrt(mmmm(iloo,1)**2.0_dp+ &
mmmm(iloo,2)**2.0_dp+mmmm(iloo,3)**2.0_dp);
end do
end if
ro=robs(:,row(coilid))

    corrfld(1)=0.0_dp
    corrfld(2)=0.0_dp
    corrfld(3)=0.0_dp

boldk3=0.0_dp
u0=(ro(1)-q(1,1))*uhat(1)+(ro(2)-q(1,2))*uhat(2)+(ro(3)-q(1,3))*uhat(3)
v0=(ro(1)-q(1,1))*vhat(1)+(ro(2)-q(1,2))*vhat(2)+(ro(3)-q(1,3))*vhat(3)
w0=(ro(1)-q(1,1))*nhat(1)+(ro(2)-q(1,2))*nhat(2)+(ro(3)-q(1,3))*nhat(3)
s_m(1)=-((l(3)-u0)*(l(3)-u3)+v0*v3)/l(1)
s_m(2)=-(u3*(u3-u0)+v3*(v3-v0))/l(2)
s_m(3)=-u0
s_p=s_m+l
t0(1)=(v0*(u3-l(3))+v3*(l(3)-u0))/l(1)
t0(2)=(u0*v3-v0*u3)/l(2)
t0(3)=v0
R0(1)=sqrt(t0(1)**2.0_dp+w0**2.0_dp)
R0(2)=sqrt(t0(2)**2.0_dp+w0**2.0_dp)
R0(3)=sqrt(t0(3)**2.0_dp+w0**2.0_dp)
Rp(1)=sqrt((ro(1)-q(3,1))**2.0_dp+(ro(2)-q(3,2))**2.0_dp+(ro(3)-q(3,3))**2.0_dp)
Rp(2)=sqrt((ro(1)-q(1,1))**2.0_dp+(ro(2)-q(1,2))**2.0_dp+(ro(3)-q(1,3))**2.0_dp)
Rp(3)=sqrt((ro(1)-q(2,1))**2.0_dp+(ro(2)-q(2,2))**2.0_dp+(ro(3)-q(2,3))**2.0_dp)
Rm(1)=Rp(3)
Rm(2)=Rp(1)
Rm(3)=Rp(2)
I(1)=(log((Rp(1)+s_p(1))/(Rm(1)+s_m(1))))
I(2)=(log((Rp(2)+s_p(2))/(Rm(2)+s_m(2))))
I(3)=(log((Rp(3)+s_p(3))/(Rm(3)+s_m(3))))

if (I(1)/=I(1)) then
I(1)=0.0_dp;
else if (I(1)-1.0_dp==I(1)) then
I(1)=0.0_dp
end if

if (I(2)/=I(2)) then
I(2)=0.0_dp;
else if (I(2)-1.0_dp==I(2)) then
I(2)=0.0_dp
end if
if (I(3)/=I(3)) then
I(3)=0.0_dp
else if (I(3)-1.0_dp==I(3)) then
I(3)=0.0_dp
end if

boldI(1)=mmmm(1,1)*I(1)+mmmm(2,1)*I(2)+mmmm(3,1)*I(3)
boldI(2)=mmmm(1,2)*I(1)+mmmm(2,2)*I(2)+mmmm(3,2)*I(3)
boldI(3)=mmmm(1,3)*I(1)+mmmm(2,3)*I(2)+mmmm(3,3)*I(3)

if (abs(w0)<=1d-12) then
do iloo=1,3
I(iloo)=0.0_dp

if ((abs(Rp(iloo)+s_p(iloo))>=1d-12).AND.(abs(Rm(iloo)+s_m(iloo))>=1d-12)) then
I(iloo)=(log((Rp(iloo)+s_p(iloo))/(Rm(iloo)+s_m(iloo))))
end if
end do


if (I(1)/=I(1)) then
I(1)=0.0_dp;
else if (I(1)-1.0_dp==I(1)) then
I(1)=0.0_dp
end if

if (I(2)/=I(2)) then
I(2)=0.0_dp;
else if (I(2)-1.0_dp==I(2)) then
I(2)=0.0_dp
end if
if (I(3)/=I(3)) then
I(3)=0.0_dp
else if (I(3)-1.0_dp==I(3)) then
I(3)=0.0_dp
end if
boldI(1)=mmmm(1,1)*I(1)+mmmm(2,1)*I(2)+mmmm(3,1)*I(3)
boldI(2)=mmmm(1,2)*I(1)+mmmm(2,2)*I(2)+mmmm(3,2)*I(3)
boldI(3)=mmmm(1,3)*I(1)+mmmm(2,3)*I(2)+mmmm(3,3)*I(3)
K1(1)=0.0_dp
K1(2)=(t0(1)*I(1)+t0(2)*I(2)+t0(3)*I(3))

else
beta(1)=atan((t0(1)*s_p(1))/(R0(1)**2.0_dp+abs(w0)*Rp(1)))-atan((t0(1)*s_m(1))/(R0(1)**2.0_dp+abs(w0)*Rm(1)))
beta(2)=atan((t0(2)*s_p(2))/(R0(2)**2.0_dp+abs(w0)*Rp(2)))-atan((t0(2)*s_m(2))/(R0(2)**2.0_dp+abs(w0)*Rm(2)))
beta(3)=atan((t0(3)*s_p(3))/(R0(3)**2.0_dp+abs(w0)*Rp(3)))-atan((t0(3)*s_m(3))/(R0(3)**2.0_dp+abs(w0)*Rm(3)))
             K1(1)=1/abs(w0)*(beta(1)+beta(2)+beta(3))
K1(2)=-(w0**2.0_dp)*K1(1)+(t0(1)*I(1)+t0(2)*I(2)+t0(3)*I(3))
end if

boldk3=-w0*nhat*K1(1)-boldI

Efldcorr(1,row(coilid))=Efldcorr(1,row(coilid))+(-boldk3(1))*xval(col(coilid))
Efldcorr(2,row(coilid))=Efldcorr(2,row(coilid))+(-boldk3(2))*xval(col(coilid))
Efldcorr(3,row(coilid))=Efldcorr(3,row(coilid))+(-boldk3(3))*xval(col(coilid))
end do
end subroutine




end module nearfielding
