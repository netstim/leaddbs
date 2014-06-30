function [coords,goodz]=ea_sample_cuboid(trajectory,trajvector,options)

switch options.modality
    case 1 % MR
        if exist([options.root,options.prefs.patientdir,filesep,options.prefs.cornii],'file')
        niifn=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii];
        else
        niifn=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii];
        end
    case 2 % CT
        niifn=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii];
end




V=spm_vol(niifn);

top_vx=[trajectory(1:2,:),ones(2,1)]';
top_mm=V.mat*top_vx;
top_vx=top_vx(1:3,:)';
top_mm=top_mm(1:3,:)';

% % map trajectory to voxels:
% trajectory_mm=[trajectory,ones(length(trajectory),1)];
%     trajectory_vx=V.mat\trajectory_mm';
%     trajectory_vx=trajectory_vx(1:3,:)';
%     trajectory_mm=trajectory_mm(:,1:3);

   
    
dvox=pdist(top_vx);
dmm=pdist(top_mm);
mm2vx=dmm/dvox; % -> 1 mm equals mm2vx voxels.
%[trajectory,trajectory]=ea_map_coords(trajectory',traniifn);
%trajectory=trajectory';
%trajectory=trajectory';
reldist=options.elspec.eldist/mm2vx; % real measured distance between electrodes in voxels.


% interpolate to include all z-heights:

% xversatz=mean(diff(trajectory(1:end,1))); %wmean(diff(centerline(1:end,1)),gaussweights,1);
% yversatz=mean(diff(trajectory(1:end,2)));%wmean(diff(centerline(1:end,2)),gaussweights,1);
% zversatz=mean(diff(trajectory(1:end,3)));%wmean(diff(centerline(1:end,2)),gaussweights,1);
% 
% trajvector2=[xversatz,yversatz,zversatz];

ntrajvector=trajvector/norm(trajvector);

trajvector=trajvector*((reldist/10)/(norm(trajvector))); % normed trajvector.
if trajvector(3)<0
    trajvector=trajvector*-1; % now going from dorsal to ventral.
    ntrajvector=ntrajvector*-1;
end


% now set trajvector to a size that 10 of it will be the size of eldist.





startpt=trajectory(end,:);            % 3d point of starting line (5 vox below coord 1).
%endpt=coords_vox(1,:);            % 3d point of end line (5 vox above coord 4).


orth=null(trajvector);
orthx=orth(:,1)'*((reldist/10)/norm(orth(:,1))); % vector going perpendicular to trajvector by 0.5 mm in x dir.
orthy=orth(:,2)'*((reldist/10)/norm(orth(:,2))); % vector going perpendicular to trajvector by 0.5 mm in y dir.




xdim=15;
ydim=15;
zdim=150; % will be sum up to 5 times reldist (three between contacts and two at borders).


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
                % plot3(coord2extract(1),coord2extract(2),coord2extract(3),'.');
            end
        end
    end
    
%                     imat(xx+xdim+1,yy+ydim+1,zz,(tracor))=spm_sample_vol(V,coord2extract(1),coord2extract(2),coord2extract(3),3);
    imat(sub2ind(size(imat),coord2write(:,1),coord2write(:,2),coord2write(:,3)))=spm_sample_vol(V,coord2extract(:,1),coord2extract(:,2),coord2extract(:,3),3);

    





%imat=ea_gencontrastimage(imat,options.zheights);

%cnii=make_nii(imat);
%save_nii(cnii,'trajectory_run.nii');

switch options.modality
    case 1
        headtemp=load_nii(fullfile(options.earoot,'templates','electrode_contacts','mr','template.nii'));
    case 2
        headtemp=load_nii(fullfile(options.earoot,'templates','electrode_contacts','ct','template.nii'));
end

htemp=headtemp.img;
temp=nanmean(htemp(:,:,:,:),4);
cimat=squeeze(imat(:,:,:));




%% calculate cross-correlation series for all x and y values within the
%% cuboid volume.

disp('Calculating cross-correlation series for all x and y values within the cuboid volume.');

    
    
 %   cnii=make_nii(temp);
 %   save_nii(cnii,['template_run',num2str(1),'.nii']);
    
cnt=1;
corrs=zeros(size(cimat,1)*size(cimat,2),length(min(trajectory(:,3)):size(cimat,3)-size(temp,3)));
for xx=1:size(cimat,1)
    for yy=1:size(cimat,2)
        for lag=min(trajectory(:,3)):size(cimat,3)-size(temp,3)
            corrs(cnt,lag+1-min(trajectory(:,3)))=corr(squeeze(cimat(xx,yy,lag:lag+size(temp,3)-1)),squeeze(temp(xx,yy,:)),'rows','pairwise');

        end
                    cnt=cnt+1;
    end
end

corrs=nanmean(corrs);
[maxv,maxi]=max(corrs);


try
goodz=maxi+min(trajectory(:,3))+reldist; % -1 because lag-series starts with 0 and not one. +reldist because templates start ~reldist voxels before electrode point.
ea_showdis(['Maximal value was: ',num2str(maxv),'.'],options.verbose);
catch
   ea_showdis(['Probably algorithm stopped to early. No guess of electrode heights possible.'],options.verbose);
   goodz=min(trajectory(:,3));
end



for coo=1:4
    coords(coo,:)=startpt+ntrajvector*maxi       +ntrajvector*(coo*reldist);

end










