function y=matvecD(Achlin,rfmm,nhatfmm,nfmm,Anear,x)
ytemp=Achlin*(x/(4*pi));
errfmm=10^-6;
mex_id_ = 'fmmD(i double[x], i double[x], i double[xx], i double[xx], i int[x], o double[x])';
[ytemp] = BEM(mex_id_, errfmm, ytemp, rfmm, nhatfmm, nfmm, 1, nfmm, 3, nfmm, 3, nfmm, 1, nfmm);
y=Achlin'*ytemp+Anear*x;
end
