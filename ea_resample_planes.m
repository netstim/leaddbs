function imat=ea_resample_planes(V,meanfitline,sample_width,doxx,resolution)
% this function samples the 2D-planes orthogonally to the fitted trajectory
% line. The code has gotten a bit cryptic for speed improvement
% reasons.

xxvec=-sample_width:resolution:sample_width;
nn=length(meanfitline);
xxlength=length(xxvec);
imat=zeros(nn,xxlength);

meanfitline=[meanfitline;ones(1,nn)]; % add ones for voxeltransformation.
addvolume=repmat(xxvec,nn,1);

fitvolume=repmat(meanfitline,[1,1,xxlength]);
if ~doxx
fitvolume(2,:,:)=squeeze(fitvolume(2,:,:))+addvolume;
elseif doxx==1
fitvolume(1,:,:)=squeeze(fitvolume(1,:,:))+addvolume;
elseif doxx==2
    fitvolume(3,:,:)=squeeze(fitvolume(3,:,:))+addvolume;
end

fitvolume=V.mat \ reshape(fitvolume,4,xxlength*nn);
flags.interp=4;
flags.wrap=[0,0,0];
d       = [flags.interp*[1 1 1]' flags.wrap(:)];


imat(:)=spm_sample_vol(V,double(fitvolume(1,:)),double(fitvolume(2,:)),double(fitvolume(3,:)),1);




function fitvolume=submax(fitvolume,sz)
for dim=1:3
fitvolume(dim,fitvolume(dim,:)<1)=1;
fitvolume(dim,fitvolume(dim,:)>sz(dim))=sz(dim);
end