function [pot,fld]=computeEandVsecondaryFMM(sigdiffV,sigdiffE,Alin,rfmm,nhatfmm,nfmm,robs)
nobs=numel(robs)/3;
Vfmm=Alin*sigdiffV;
Efmm=Alin*sigdiffE;
errfmm=10^-6;
mex_id_ = 'fmmDEsecond(i double[x], i double[x], i double[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], o double[x], o double[xx])';
[pot, fld] = BEM(mex_id_, errfmm, Vfmm, Efmm, rfmm, nhatfmm, nfmm, robs, nobs, 1, nfmm, nfmm, 3, nfmm, 3, nfmm, 1, 3, nobs, 1, nobs, 3, nobs);
pot=pot/(4*pi);
fld=-fld/(4*pi);
end

