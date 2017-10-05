function ea_corrnifti(fis,regressor,output)
    % creates a correlation nifti file given a set of images and a
    % regressor
   
    
    for f=1:length(fis)
       n=ea_load_nii(fis{f});
       if ~exist('X','var')
           X=n.img(:);
       else
           X=[X,n.img(:)];
       end
    end

    R=corr(regressor,X','type','Spearman');

n.img(:)=R;
n.fname=output;

n.dt=[16,0];
ea_write_nii(n);
    