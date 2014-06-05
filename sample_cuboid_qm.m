function [goodx,goody,goodz]=sample_cuboid_qm(coords,fitline,patientname,options,side)

traniifn=[options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii'];
corniifn=[options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt.nii'];


[fitline,fitline_vox]=map_coords(fitline',traniifn);
fitline=fitline';
fitline_vox=fitline_vox';

xdim=5;
ydim=5;
zdim=size(fitline_vox,1);


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
    
    xversatz=mean(diff(fitline_vox(1:end,1))); %wmean(diff(centerline(1:end,1)),gaussweights,1);
    yversatz=mean(diff(fitline_vox(1:end,2)));%wmean(diff(centerline(1:end,2)),gaussweights,1);
    zversatz=mean(diff(fitline_vox(1:end,3)));%wmean(diff(centerline(1:end,2)),gaussweights,1);
    
    trajvector=[xversatz,yversatz,zversatz];
    
    orth=null(trajvector);
    orthx=orth(:,1);
    orthy=orth(:,2);
    
    
    
    
    for zz=1:zdim
        for xx=-xdim:xdim
            for yy=-ydim:ydim
                
                
                    coord2extract=[fitline_vox(zz,1)+orthx(1)*xx+orthy(1)*yy; ...
                        fitline_vox(zz,2)+orthx(2)*xx+orthy(2)*yy; ...
                        fitline_vox(zz,3)+orthx(3)*xx+orthy(3)*yy];
                    
                    
                    imat(xx+xdim+1,yy+ydim+1,fitline_vox(zz,3),(tracor+1))=spm_sample_vol(V,coord2extract(1),coord2extract(2),coord2extract(3),3);
         % plot3(coord2extract(1),coord2extract(2),coord2extract(3),'.');      
            end
        end
    end
    
    
    
end

[coords,coords_vox]=map_coords(coords',traniifn);
coords=coords'; coords_vox=coords_vox';

imat_crop=imat(:,:,round(coords_vox(1,3))-5:round(coords_vox(4,3))+5,:);

cnii=make_nii(imat_crop);


if options.eldist==3
    save_nii(cnii,[pwd,'/qm/3mm/',patientname,side,'.nii']);
elseif options.eldist==2
    save_nii(cnii,[pwd,'/qm/2mm/',patientname,side,'.nii']);
end


