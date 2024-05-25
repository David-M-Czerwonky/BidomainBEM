function y=matvecDoublelayer(NttoN,Cdiff,Alin,rfmm,nhatfmm,nfmm,Anear,x,precon)
x(end+1)=0;
%Doublelayer=ID/2-NttoN'*D*Cdiff*NttoN;
y=Anear*x;
ytemp=Alin*Cdiff*NttoN*(x/(4*pi));
errfmm=10^-6;
mex_id_ = 'fmmD(i double[x], i double[x], i double[xx], i double[xx], i int[x], o double[x])';
[ytemp] = BEM(mex_id_, errfmm, ytemp, rfmm, nhatfmm, nfmm, 1, nfmm, 3, nfmm, 3, nfmm, 1, nfmm);
y=y-NttoN'*Alin'*ytemp;
y=precon(1:end-1).*y(1:end-1)+precon(end)*y(end)/sqrt(numel(x));
end

