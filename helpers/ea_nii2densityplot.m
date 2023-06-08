function h=ea_nii2densityplot(niifname,color,subsample,resultfig)
if ~exist('resultfig','var')
    resultfig=figure;
end


% if subsample
%     t=ea_getleadtempdir;
%     uid=ea_generate_uuid;
%     copyfile(niifname,[t,uid,'.nii']);
%     ea_reslice_nii([t,uid,'.nii'],[t,uid,'.nii'],[subsample,subsample,subsample]);
%     niifname=[t,uid,'.nii'];
% end

nii=ea_load_nii(niifname);



[xx,yy,zz]=ind2sub(size(nii.img),1:numel(nii.img));
XYZ=nii.mat*[xx;yy;zz;ones(1,length(xx))];
XYZ=XYZ(1:3,:)';
voxsz=mean(nii.voxsize);
intensities=ea_minmax(nii.img(:));
ea_dispercent(0,'Plotting spheres')
set(0,"CurrentFigure",resultfig)

samples=1:subsample:length(XYZ);
samples=round(samples+(randn(1,size(samples,2)).*subsample));
if samples(1)<1
    samples(1)=1;
end
if samples(end)>length(XYZ)
    samples(end)=length(XYZ);
end
samples(isnan(intensities(samples)))=[];
cnt=1;
samplelen=length(samples);
for pt=samples
    usealpha=exp(4*intensities(pt)-4)./2;
        h(cnt)=ea_plotsphere(XYZ(pt,:),intensities(pt)*voxsz*0.1*subsample,color,'none',usealpha);
        cnt=cnt+1;
        %drawnow
        ea_dispercent(cnt/samplelen);
end


[h.DiffuseStrength]=deal(0.3);
[h.AmbientStrength]=deal(0.6);
[h.SpecularStrength]=deal(0);
drawnow
ea_dispercent(1,'end');
