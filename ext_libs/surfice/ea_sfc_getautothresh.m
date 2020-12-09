function threshs=ea_sfc_getautothresh(fis)


    for fi=1:length(fis)
        nii=ea_load_nii(fis{fi});
        nii.img(nii.img==0)=nan;
        thisstd=ea_nanstd(nii.img(:));
        if min(nii.img(:))<0 % has negative values, calculate 4 vals
            threshs(fi,1)=0.5*thisstd;
            threshs(fi,2)=2.7*thisstd;
            threshs(fi,3)=-0.5*thisstd;
            threshs(fi,4)=-2.7*thisstd;
        else % only positive vals
            threshs(fi,1)=0.04*thisstd;
            threshs(fi,2)=1.9*thisstd;
            threshs(fi,3)=nan;
            threshs(fi,4)=nan;
        end
        
    end
    if ~any(nii.img(:)>0)
       threshs(fi,1:2)=nan; 
    end
