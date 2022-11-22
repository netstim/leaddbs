function [goodx,goody,goodz]=ea_export_templates(coords,trajectory,patientname,options,side)
% this function can export templates based on manual measurements of the
% electrode tips and is originaly based on ea_sample_cuboid.
% uses map_coords authored by Ged Ridgway

switch options.subj.postopModality
    case 'MRI'
        if isfield(options.subj.coreg.anat.postop,'cor_MRI') && isfile(options.subj.coreg.anat.postop.cor_MRI)
            niifn=options.subj.coreg.anat.postop.cor_MRI;
        else
            niifn=isfile(options.subj.coreg.anat.postop.ax_MRI);
        end
    case 'CT'
        niifn=options.subj.coreg.anat.postop.CT;
end

[trajectory,trajectory_vox]=ea_map_coords(trajectory',niifn);
trajectory=trajectory';
trajectory_vox=trajectory_vox';

% interpolate to include all z-heights:
[coords,coords_vox]=ea_map_coords(coords',niifn);
coords=coords'; coords_vox=coords_vox';
reldist=ea_pdist(coords_vox(2:3,:)); % real measured distance between electrodes in voxels.

xversatz=mean(diff(trajectory_vox(1:end,1))); %wmean(diff(centerline(1:end,1)),gaussweights,1);
yversatz=mean(diff(trajectory_vox(1:end,2)));%wmean(diff(centerline(1:end,2)),gaussweights,1);
zversatz=mean(diff(trajectory_vox(1:end,3)));%wmean(diff(centerline(1:end,2)),gaussweights,1);

trajvector=[xversatz,yversatz,zversatz];
trajvector=trajvector*((reldist/10)/norm(trajvector)); % traversing dir_vector of 0.5 mm in distance along trajectory.
if trajvector(3)<0
    trajvector=trajvector*-1; % now going from dorsal to ventral.
end

startpt=coords_vox(1,:)-10*trajvector;            % 3d point of starting line (5 vox below coord 1).

orth=null(trajvector);
orthx=orth(:,1)'*((reldist/10)/norm(orth(:,1))); % vector going perpendicular to trajvector by 0.5 mm in x dir.
orthy=orth(:,2)'*((reldist/10)/norm(orth(:,2))); % vector going perpendicular to trajvector by 0.5 mm in y dir.

xdim=15;
ydim=15;
zdim=50; % will be sum up to 5 times reldist (three between contacts and two at borders).

imat=nan(2*ydim+1,2*xdim+1,zdim);
V=spm_vol(niifn);

cnt=1;
coord2write=zeros(length(1:zdim)* length(-xdim:xdim)*length(-ydim:ydim),3);
coord2extract=zeros(length(1:zdim)* length(-xdim:xdim)*length(-ydim:ydim),3);

for zz=1:zdim
    for xx=-xdim:xdim
        for yy=-ydim:ydim

            pt=startpt+zz*trajvector;
            coord2extract(cnt,:)=[pt(1)+orthx(1)*xx+orthy(1)*yy; ...
                pt(2)+orthx(2)*xx+orthy(2)*yy; ...
                pt(3)+orthx(3)*xx+orthy(3)*yy]';
            coord2write(cnt,:)=[xx+xdim+1;yy+ydim+1;zz]';
            cnt=cnt+1;
        end
    end
end

imat(sub2ind(size(imat),coord2write(:,1),coord2write(:,2),coord2write(:,3)))=spm_sample_vol(V,coord2extract(:,1),coord2extract(:,2),coord2extract(:,3),3);

% save templates.
switch options.subj.postopModality
    case 'MRI'
        mrstr='mr';
    case 'CT'
        mrstr='ct';
end

if exist([ea_getearoot,'templates',filesep,'electrode_contacts',filesep,mrstr,filesep,'template.nii'],'file')
    template=load_nii([ea_getearoot,'templates',filesep,'electrode_contacts',filesep,mrstr,filesep,'template.nii']);
    nutimg=zeros(size(template.img,1),size(template.img,2),size(template.img,3),size(template.img,4)+1);
    nutimg(:,:,:,1:end-1)=template.img;
else
    nutimg=zeros(2*xdim+1,2*ydim+1,zdim,1);
end

nutimg(:,:,:,end)=imat;
cnii=make_nii(nutimg);
save_nii(cnii,[ea_getearoot,'templates',filesep,'electrode_contacts',filesep,mrstr,filesep,'template.nii']);
