function [goodx,goody,goodz]=sample_cuboid(diams,fitline,trajvector,patientname,zdist,options)

traniifn=[options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii'];
corniifn=[options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt.nii'];


xdim=5;
ydim=5;
zdim=max(fitline(:,3));


imat=zeros(2*ydim+1,2*xdim+1,zdim,2);


for tracor=0:1
    
    if tracor==1
        nii=load_nii(corniifn); %'_cor_brain_A3_final_wo_opt.nii']);
        V=spm_vol(corniifn);
    elseif tracor==0
        nii=load_nii(traniifn); %'_cor_brain_A3_final_wo_opt.nii']);
        V=spm_vol(traniifn);
    end
    nii.img=double(nii.img);
    
    img=nii.img;
    
    
    
    
    
    orth=null(trajvector);
    orthx=orth(:,1);
    orthy=orth(:,2);
    
    
    
    
    
    for zz=1:zdim
        for xx=-xdim:xdim
            for yy=-ydim:ydim
                
                try
                    coord2extract=[fitline(zz,1)+orthx(1)*xx+orthy(1)*yy; ...
                        fitline(zz,2)+orthx(2)*xx+orthy(2)*yy; ...
                        fitline(zz,3)+orthx(3)*xx+orthy(3)*yy];
                    
                    
                    imat(xx+xdim+1,yy+ydim+1,fitline(zz,3),(tracor+1))=spm_sample_vol(V,coord2extract(1),coord2extract(2),coord2extract(3),3);
                end
            end
        end
    end
    
    
    
end



imat=gencontrastimage(imat,options.zheights);





bestval=inf; % high value at start.
onsetvolume=nan(size(imat,1),size(imat,2),length(1:(1/options.zresolution):length(imat)));
for xx=1:size(imat,1)
    
    for yy=1:size(imat,2)
        
        
        % build new trace built up from diameters of image and imat line.
        
        conv_diams=squeeze(imat(xx,yy,:));
        %conv_diams(1:length(diams))=conv_diams(1:length(diams)).*diams'; % multiply with old diams..
        
        %conv_diams(length(diams):end)=conv_diams(length(diams):end).*mean(diams(diams>0));
        
        %diams=squeeze(imat(xx,yy,:));
        
        
        
        onsetvolume(xx,yy,:)=match_humps(conv_diams,zdist,options);
        
    end
end


% this new approach will overwrite goodz based on
nii=make_nii(onsetvolume);
save_nii(nii,'onsetvolume.nii');
spm_smooth('onsetvolume.nii','onsetvolume.nii',[4 4 4]);
nii=load_nii('onsetvolume.nii');

nii.img(:,:,1:200)=-inf; % cut out these z-regions
nii.img(:,:,500:end)=-inf;

hdfitline=genhdfitline(fitline,options);

% delete all entries that are not covered by hdfitline.
tmpnii=inf(size(nii.img,1),size(nii.img,2),size(nii.img,3));
tmpnii=tmpnii*-1;

diamidx=1:length(conv_diams);
diamidx=interp1q([1:length(diamidx)]',diamidx',[1:1/options.zresolution:length(diamidx)]');
minheight=find(diamidx==min(hdfitline(:,3)));
maxheight=find(diamidx==max(hdfitline(:,3)));

tmpnii(:,:,[minheight:maxheight])=0;
nii.img=nii.img+tmpnii;



if options.priorstrength;
    genpriorvolume(options);
    priornii=load_nii('prior.nii');
    nii.img=((1-options.priorstrength)*nii.img).*(options.priorstrength*priornii.img);
end





[goodx,goody,goodz]=ind2sub(size(nii.img),find(nii.img(:)==max(nii.img(:))));
goodz=goodz(1); % for now if the rare case happens of 2 minima, just guess and take the first one.
goodx=goodx-size(nii.img,1);
goody=goody-size(nii.img,2);

showdis(['Maximal value was: ',num2str(max(nii.img(:))),'.'],options.verbose);


%% re-downsample to voxel space


goodz=diamidx(goodz); % resample to downsampled "real" voxel values.



% if options.targetknown
%     switch lower(options.target)
%         case 'stn'
%             canonicalz=30;
%         case 'gpi'
%             canonicalz=40;
%         case 'vim'
%             canonicalz=50;
%     end
% [dummy,goodzs]=findpeaks(allonsetmap.*-1);
%
% diffs=diff([goodzs;repmat(canonicalz,1,length(goodzs))]);
% goodz=goodzs(abs(diffs)==min(abs(diffs)));
% -> gives out the local minimum in allonsetmap closest to the prior
% canonicalz.
% end





%% visualize
if options.verbose>1;
    discor=figure;
    
    imagesc(mean(imat,3));
    
    
    colormap gray
    
    
    
    
    if options.verbose>2;
        close(discor);
    end
end




% function vox=tra2cor(vox,traniifn,corniifn,tracor)
%
% if tracor
%     keyboard
%     vox=map_coords(vox,traniifn); % -> to world-coordinates
%     [dummy,vox]=map_coords(vox,corniifn); % -> to cor-voxelspace
% end
%
% vox=round(vox);




