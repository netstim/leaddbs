function [coords,goodz]=ea_reconstruct_coords(trajectory,trajvector,options)

[cimat,reldist]=ea_sample_cuboid(trajectory,options);

vizz=0;

switch options.subj.postopModality
    case 'MRI'
        headtemp=ea_load_nii(fullfile(ea_getearoot,'templates',filesep,'electrode_contacts','mr','template.nii'));
    case 'CT'
        headtemp=ea_load_nii(fullfile(ea_getearoot,'templates',filesep,'electrode_contacts','ct','template.nii'));
end

htemp=headtemp.img;
temp=ea_nanmean(htemp(:,:,:,:),4);

% calculate cross-correlation series for all x and y values within the cuboid volume.
disp('Calculating cross-correlation series for all x and y values within the cuboid volume.');

cnt=1;
corrs=zeros(size(cimat,1)*size(cimat,2),299);
for xx=1:size(cimat,1)
    for yy=1:size(cimat,2)
        [corrs(cnt,:),lags]=ea_xcorr(squeeze(cimat(xx,yy,:)),squeeze(temp(xx,yy,:)));
        cnt=cnt+1;
    end
end

corrs=ea_nanmean(corrs);
[maxv,maxi]=max(corrs);
maxi=lags(maxi);

if vizz
   figure
   subplot(2,1,1);
   icimat=cat(3,cimat,zeros(31,31,maxi));
   imagesc(squeeze(icimat(15,:,:)));
   colormap(gray)

   subplot(2,1,2);
   itemp=cat(3,zeros(31,31,maxi),temp);
   imagesc(squeeze(itemp(15,:,:)));
   colormap(gray)
end

try
    goodz=maxi+min(trajectory(:,3))+reldist; % -1 because lag-series starts with 0 and not one. +reldist because templates start ~reldist voxels before electrode point.
    ea_showdis(['Maximal value was: ',num2str(maxv),'.'],options.verbose);
catch
   ea_showdis(['Probably algorithm stopped to early. No guess of electrode heights possible.'],options.verbose);
   goodz=min(trajectory(:,3));
   maxi=1;
end

if isempty(goodz)
	ea_showdis(['Probably algorithm stopped to early. No guess of electrode heights possible.'],options.verbose);
    goodz=min(trajectory(:,3));
    maxi=1;
end

% error on maxi empty:
if isempty(maxi)
    ea_error('Reconstruction failed. Please adjust mask-size and entry-point parameters and re-run reconstruction.');
end

startpt=trajectory(end,:);
ntrajvector=trajvector/norm(trajvector);

for coo=1:4
    coords(coo,:)=startpt-ntrajvector*maxi       -ntrajvector*((coo-1)*reldist);
end
