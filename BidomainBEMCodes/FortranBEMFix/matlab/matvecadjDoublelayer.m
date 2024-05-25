function y=matvecadjDoublelayer(NttoN,Cdiff,Alin,rfmm,nhatfmm,nfmm,Anear,x,precon)
%Doublelayer=-ID/2+NttoN'*Cdiff*Dstar*NttoN;
y=Anear*x;
ytemp=Alin*NttoN*(x/(4*pi));
errfmm=10^-6;
mex_id_ = 'fmmDstar(i double[x], i double[x], i double[xx], i double[xx], i int[x], o double[x])';
[ytemp] = BEM(mex_id_, errfmm, ytemp, rfmm, nhatfmm, nfmm, 1, nfmm, 3, nfmm, 3, nfmm, 1, nfmm);
y=y+NttoN'*Cdiff*Alin'*ytemp;
y=precon(1:end).*y(1:end);
end
