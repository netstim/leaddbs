function imat=ea_resample_planes(Vcor,meanfitline,sample_width,doxx,resolution)
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
else
fitvolume(1,:,:)=squeeze(fitvolume(1,:,:))+addvolume;
end

fitvolume=Vcor.mat \ reshape(fitvolume,4,xxlength*nn);

try % has been loaded using nifti command
    fitvolume=submax(fitvolume,size(Vcor.dat));
    imat(:)=Vcor.dat(sub2ind(size(Vcor.dat),fitvolume(1,:),fitvolume(2,:),fitvolume(3,:)));
catch % has been loaded using spm_vol command
    imat(:)=spm_sample_vol(Vcor,double(fitvolume(1,:)),double(fitvolume(2,:)),double(fitvolume(3,:)),1);
end





function fitvolume=submax(fitvolume,sz)
for dim=1:3
fitvolume(dim,fitvolume(dim,:)<1)=1;
fitvolume(dim,fitvolume(dim,:)>sz(dim))=sz(dim);
end