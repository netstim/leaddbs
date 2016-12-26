function iscoreg=ea_hdr_iscoreg(Vcomp,Vref)
iscoreg=1;

if ~isequal(Vref.dim,Vcomp.dim)
    iscoreg=0;
    return
end
if any(abs(Vref.mat(:)-Vcomp.mat(:))>0.0001) % very small rounding errors should be ignored, thus not use isequal here.
    iscoreg=0;
    return
end