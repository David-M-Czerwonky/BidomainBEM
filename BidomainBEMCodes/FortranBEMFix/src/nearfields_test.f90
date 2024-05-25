


subroutine nearmat(p,np,t2p,nt,col,row,ncols,normal,area,Dout,Sout,Doutlin,Soutlin)
  use globalconstants
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
real(kind=dp),dimension(3):: corrfld,corrtmp,cen,l,nhat,uhat,vhat,s_m,s_p,t0,Rp
real(kind=dp),dimension(3):: Rm,I,boldI,boldk3,beta,K1,R0,f3i,bk3,boldUa,boldVa
real(kind=dp),dimension(3,3):: dir,mmmm,q,f_i,f1,f2

Sout=0.0_dp
Dout=0.0_dp
Soutlin=0.0_dp
Doutlin=0.0_dp

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

f3i(1)=s_p(1)*Rp(1)-s_m(1)*Rm(1)+R0(1)**2.0_dp*I(1)
f3i(2)=s_p(2)*Rp(2)-s_m(2)*Rm(2)+R0(2)**2.0_dp*I(2)
f3i(3)=s_p(3)*Rp(3)-s_m(3)*Rm(3)+R0(3)**2.0_dp*I(3)

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

N10=1.0_dp-u0/l(3)+(u3/l(3)-1.0_dp)*v0/v3
N20=u0/l(3)-(u3/l(3))*v0/v3
N30=v0/v3


f1(:,1)=f1(:,1)+qwt(jquad)*qpt(:,jquad)* &
((N10*bk3(1)-boldUa(1)/l(3)+(u3/l(3)-1.0_dp)*boldVa(1)/v3)*nhat(1)+ &
 (N10*bk3(2)-boldUa(2)/l(3)+(u3/l(3)-1.0_dp)*boldVa(2)/v3)*nhat(2)+ &
 (N10*bk3(3)-boldUa(3)/l(3)+(u3/l(3)-1.0_dp)*boldVa(3)/v3)*nhat(3))
f1(:,2)=f1(:,2)+qwt(jquad)*qpt(:,jquad)* &
((N20*bk3(1)+boldUa(1)/l(3)-u3/l(3)*boldVa(1)/v3)*nhat(1)+ &
 (N20*bk3(2)+boldUa(2)/l(3)-u3/l(3)*boldVa(2)/v3)*nhat(2)+ &
 (N20*bk3(3)+boldUa(3)/l(3)-u3/l(3)*boldVa(3)/v3)*nhat(3))
f1(:,3)=f1(:,3)+qwt(jquad)*qpt(:,jquad)* &
((N30*bk3(1)+boldVa(1)/v3)*nhat(1)+ &
 (N30*bk3(2)+boldVa(2)/v3)*nhat(2)+ &
 (N30*bk3(3)+boldVa(3)/v3)*nhat(3))

f2(:,1)=f2(:,1)+qwt(jquad)*qpt(:,jquad)*(N10*K1(2)-Iu/l(3)+(u3/l(3)-1.0_dp)*Iv/v3)
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
area(row(coilid))*((boldk3(1))*normal(1,row(coilid))+ &
                   (boldk3(2))*normal(2,row(coilid))+ &
                   (boldk3(3))*normal(3,row(coilid)))
!Dout2(:,coilid)=Dout2(:,coilid)+ &
!area(row(coilid))*boldk3(:)

Doutlin(coilid,1)=Doutlin(coilid,1)-area(row(coilid))*f1(1,1)
Doutlin(coilid,2)=Doutlin(coilid,2)-area(row(coilid))*f1(2,1)
Doutlin(coilid,3)=Doutlin(coilid,3)-area(row(coilid))*f1(3,1)
Doutlin(coilid,4)=Doutlin(coilid,4)-area(row(coilid))*f1(1,2)
Doutlin(coilid,5)=Doutlin(coilid,5)-area(row(coilid))*f1(2,2)
Doutlin(coilid,6)=Doutlin(coilid,6)-area(row(coilid))*f1(3,2)
Doutlin(coilid,7)=Doutlin(coilid,7)-area(row(coilid))*f1(1,3)
Doutlin(coilid,8)=Doutlin(coilid,8)-area(row(coilid))*f1(2,3)
Doutlin(coilid,9)=Doutlin(coilid,9)-area(row(coilid))*f1(3,3)

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
subroutine precorr(triafl,nt,nquadin,col,row,ncols,normal,Sout,Dout)
  use globalconstants
  implicit none
integer(kind=spint),intent(in) :: nt,nquadin,ncols
integer(kind=spint),intent(in),dimension(ncols) :: col,row
real(kind=dp),intent(in),dimension(3,nquadin*nt) :: triafl
real(kind=dp),intent(in),dimension(3,nt):: normal
real(kind=dp),intent(out),dimension(ncols,nquadin,nquadin):: Sout,Dout!,Dstarout

!internal variables
real(kind=dp),dimension(3,nquadin) :: rs,ro
integer(kind=spint) :: coilid,jquad,iquad,xpre,indtemp,jndtemp
real(kind=dp):: corrpot,radtmp,radtmp1,corrD0,corrDstar0
real(kind=dp),dimension(3):: corrtmp

xpre=0


do coilid=1,ncols
if (xpre/=col(coilid)) then
xpre=col(coilid)
!get source points
do iquad=1,nquadin
rs(:,iquad)=triafl(:,iquad+(col(coilid)-1)*nquadin)
end do
end if
!get observation points
do iquad=1,nquadin
ro(:,iquad)=triafl(:,iquad+(row(coilid)-1)*nquadin)
end do
!compute correction
  corrpot=0.0_dp
  corrD0=0.0_dp
  corrDstar0=0.0_dp
do jquad=1,nquadin
do iquad=1,nquadin
corrtmp(:)=ro(:,jquad)-rs(:,iquad)
radtmp=sqrt(corrtmp(1)**2+corrtmp(2)**2+corrtmp(3)**2)
if (radtmp>=1d-14) then
  indtemp=iquad+nquadin*(col(coilid)-1)
  jndtemp=jquad+nquadin*(row(coilid)-1)
  radtmp1=1.0_dp/radtmp
  radtmp=1.0_dp/radtmp**3
  !corrpot=corrpot-radtmp1
  Sout(coilid,jquad,iquad)=radtmp1
!  radtmp1=(corrtmp(1)*normal(1,row(coilid))+ &
!            corrtmp(2)*normal(2,row(coilid))+ &
!            corrtmp(3)*normal(3,row(coilid)))*radtmp
!corrDstar0=corrDstar0-radtmp1
!Dstarout(coilid,jquad,iquad)=-radtmp1
radtmp1=(corrtmp(1)*normal(1,col(coilid))+ &
          corrtmp(2)*normal(2,col(coilid))+ &
          corrtmp(3)*normal(3,col(coilid)))*radtmp
!corrD0=corrD0-radtmp1
Dout(coilid,jquad,iquad)=radtmp1
else
  Sout(coilid,jquad,iquad)=0.0_dp
Dout(coilid,jquad,iquad)=0.0_dp
end if
end do
end do
!save correction



end do


end subroutine


subroutine fmmS(errfmm,chspace,rfmm,nfmm,pot)
  use globalconstants
implicit none
real(kind=dp),intent(in) :: errfmm
real(kind=dp),dimension(3,nfmm),intent(in) :: rfmm
real(kind=dp),dimension(nfmm),intent(in) :: chspace
real(kind=dp),dimension(nfmm),intent(out) :: pot
integer(kind=spint),intent(in) :: nfmm
integer(kind=spint) :: flag

          call lfmm3d_s_c_p(errfmm,nfmm,rfmm,chspace,pot)
!if flag is 4/8 not enough memory flag=0 means its correct
end subroutine fmmS


subroutine fmmDstar(errfmm,chspace,rfmm,nhatfmm,nfmm,pot)
  use globalconstants
implicit none
real(kind=dp),intent(in) :: errfmm
real(kind=dp),dimension(3,nfmm),intent(in) :: rfmm,nhatfmm
real(kind=dp),dimension(nfmm),intent(in) :: chspace
real(kind=dp),dimension(nfmm),intent(out) :: pot
integer(kind=spint),intent(in) :: nfmm
real(kind=dp),dimension(3,nfmm) :: fld

          call lfmm3d_s_c_g(errfmm,nfmm,rfmm,chspace,pot,fld)
pot=nhatfmm(1,:)*fld(1,:)+nhatfmm(2,:)*fld(2,:)+nhatfmm(3,:)*fld(3,:)

!if flag is 4/8 not enough memory flag=0 means its correct
end subroutine fmmDstar


subroutine fmmD(errfmm,chspace,rfmm,nhatfmm,nfmm,pot)
  use globalconstants
implicit none
real(kind=dp),intent(in) :: errfmm
real(kind=dp),dimension(3,nfmm),intent(in) :: rfmm,nhatfmm
real(kind=dp),dimension(nfmm),intent(in) :: chspace
real(kind=dp),dimension(nfmm),intent(out) :: pot
integer(kind=spint),intent(in) :: nfmm
real(kind=dp),dimension(3,nfmm) :: fld
fld(1,:)=nhatfmm(1,:)*chspace
fld(2,:)=nhatfmm(2,:)*chspace
fld(3,:)=nhatfmm(3,:)*chspace
          call lfmm3d_s_d_p(errfmm,nfmm,rfmm,fld,pot)


!if flag is 4/8 not enough memory flag=0 means its correct
end subroutine fmmD



subroutine generategroupct(t2p,nt,p,np,del,ncol)
  use globalconstants
  use neargrouping
implicit none
integer(kind=spint),intent(in) :: nt,np
integer(kind=spint),intent(in),dimension(3,nt) :: t2p
real(kind=dp),intent(in),dimension(3,np) :: p
real(kind=dp), intent(in) :: del
integer(kind=spint),intent(out) :: ncol
!internals
real(kind=dp),dimension(:,:),allocatable :: rs
integer(kind=spint),dimension(:),allocatable :: coltt,rowtt,grs
integer(kind=spint) :: ct,ct2,iloo,jloo,kloo,ns
allocate(rs(3,nt))
allocate(grs(nt))
!create groups
do iloo=1,nt
rs(:,iloo)=p(:,t2p(1,iloo))+p(:,t2p(2,iloo))+p(:,t2p(3,iloo))
end do
rs=rs/3.0_dp
ns=nt
grs=(/ (iloo, iloo = 1, ns) /)
ncol=3*nt
allocate(coltt(ncol))
allocate(rowtt(ncol))
ct=0
ct2=0
call countneighbors(t2p,nt,coltt,rowtt,ct,ct2)
!figure out how many nearfield
call groupcreatorMLct(rs,grs,ns,rs,grs,ns,del,ct)
	  ncol=ct+ct2

deallocate(coltt)
deallocate(rowtt)

end subroutine




subroutine generategroupmatlab(t2p,nt,p,np,del,colt,rowt,ncol)
use globalconstants
use neargrouping
implicit none
integer(kind=spint),intent(in) :: nt,np
integer(kind=spint),intent(in),dimension(3,nt) :: t2p
integer(kind=spint),dimension(ncol),intent(out) :: colt,rowt
real(kind=dp),intent(in),dimension(3,np) :: p
real(kind=dp), intent(in) :: del
integer(kind=spint),intent(in) :: ncol

!internals
real(kind=dp),dimension(:,:),allocatable :: rs
integer(kind=spint),dimension(:),allocatable :: coltt,rowtt,grs
integer(kind=spint) :: ct,ct2,iloo,jloo,kloo,ns
allocate(rs(3,nt))
allocate(grs(nt))
do iloo=1,nt
rs(:,iloo)=p(:,t2p(1,iloo))+p(:,t2p(2,iloo))+p(:,t2p(3,iloo))
end do
rs=rs/3.0_dp
ns=nt
grs=(/ (iloo, iloo = 1, ns) /)
allocate(coltt(3*nt))
allocate(rowtt(3*nt))
ct2=0
ct=0
call countneighbors(t2p,nt,coltt,rowtt,ct,ct2)
ct=0
call writeneighbors(t2p,nt,coltt,rowtt,colt,rowt,ncol,ct,ct2)
!write nearfield
call groupcreatorML(rs,grs,ns,rs,grs,ns,del,ct,colt,rowt,ncol)
call QsortC(colt,rowt)
deallocate(rowtt)
deallocate(coltt)
allocate(rowtt(ncol))
allocate(coltt(ncol))
coltt=0
rowtt=0
ct=1
coltt(1)=colt(1)
rowtt(1)=rowt(1)
do iloo=2,ncol
if ((colt(iloo) .ne. colt(iloo-1)) .or. (rowt(iloo) .ne. rowt(iloo-1))) then
ct=ct+1
coltt(ct)=colt(iloo)
rowtt(ct)=rowt(iloo)
end if
end do

colt=coltt
rowt=rowtt

deallocate(rowtt)
deallocate(coltt)


deallocate(rs)
deallocate(grs)


end subroutine


subroutine computeEprimary(errfmm,rs,js,ns,robs,nobs,Eprimary)
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: spint=kind(8)
  real(kind=dp),parameter :: pi=3.141592653589793_dp,epsthresh=1d-12
  real(kind=dp),parameter ::eps0=8.854187820000000d-12
  real(kind=dp) ::  mu0over4pi=-1.000000000000000d-7
  real(kind=dp),intent(in) :: errfmm
real(kind=dp),dimension(3,ns),intent(in) :: rs
real(kind=dp),dimension(3,nobs),intent(in) :: robs
real(kind=dp),dimension(3,ns),intent(in) :: js
real(kind=dp),dimension(3,nobs),intent(out) :: Eprimary
integer(kind=spint),intent(in) :: ns,nobs

integer(kind=spint) :: i,j,ier

call lfmm3d_t_c_p_vec(3,errfmm,ns,rs,js,nobs,robs,Eprimary)
Eprimary=mu0over4pi*Eprimary


end subroutine



subroutine fmmDEsecond(errfmm,voltage,Eprimdotn,rfmm,nhatfmm,nfmm,robs,nobs,pottarg,gradtarg)
  use globalconstants
implicit none
integer(kind=spint),intent(in) :: nfmm,nobs
real(kind=dp),intent(in) :: errfmm
real(kind=dp),dimension(3,nobs),intent(in) :: robs
real(kind=dp),dimension(3,nfmm),intent(in) :: rfmm,nhatfmm
real(kind=dp),dimension(nfmm),intent(in) :: Eprimdotn,voltage
real(kind=dp),dimension(nobs),intent(out) :: pottarg
real(kind=dp),dimension(3,nobs),intent(out) :: gradtarg
real(kind=dp),dimension(3,nfmm) :: fld,dipvec

dipvec(1,:)=nhatfmm(1,:)*voltage
dipvec(2,:)=nhatfmm(2,:)*voltage
dipvec(3,:)=nhatfmm(3,:)*voltage
call  lfmm3d_t_cd_g(errfmm,nfmm,rfmm,Eprimdotn,dipvec,nobs,robs,pottarg,gradtarg)
!if flag is 4/8 not enough memory flag=0 means its correct
end subroutine
