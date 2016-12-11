function ea_resliceanat(options)

directory=[options.root,options.patientname,filesep];

V=spm_vol([directory,options.prefs.prenii_unnormalized]);
probepnt=ones(4,1);
vsz=zeros(3,1);
for dim=1:3
    probepnt1=probepnt;
    probepnt1(dim)=probepnt1(dim)+1;
    pnt=V.mat*probepnt;
    pnt1=V.mat*probepnt1;
    vsz(dim)=pdist([pnt(1:3),pnt1(1:3)]');
end

if any(vsz>1)
    nvsz=vsz;
    nvsz(vsz>1)=1;
    ea_reslice_nii([directory,options.prefs.prenii_unnormalized],[directory,options.prefs.prenii_unnormalized],nvsz);
end
