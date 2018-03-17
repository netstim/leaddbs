function [X,n]=ea_genX(fis,regressor,output,mask,sk,ctype)

    if ~exist('mask','var')
       mask=nan; 
    else
       if isempty(mask)
           mask=nan;
       end
    end
    
    % check wheter want to normalize or smooth data
    if exist('sk','var')
        switch sk
            case 'k'
                dok=1; dos=0;
            case 's'
                dos=1; dok=0;
            case 'sk'
                dos=1; dok=1;
            otherwise
                dos=0; dok=0;
        end
    else
        dos=0;
        dok=0;
    end
   
    tmpd=ea_getleadtempdir;
    for f=1:length(fis)
        if dok || dos
            guid=ea_generate_guid;
            [pth,fn,ext]=fileparts(fis{f});
            if strcmp(ext,'.gz')
                addgz='.gz';
            else
                addgz='';
            end
            copyfile(fis{f},[tmpd,guid,'.nii',addgz]);
            if ~isempty(addgz)
                gunzip([tmpd,guid,'.nii',addgz]);
            end
            nii=ea_load_nii([tmpd,guid,'.nii']);
            if isnan(mask)
                mask=1:numel(nii.img);
                warning('If using k option should apply a mask');
            end
            nii.img(~mask)=nan;
            if dok
                nii.img(mask)=ea_normal(nii.img(mask),1,0,' ',0,1,'TRUE');
            end
            ea_write_nii(nii);
            if dos
                matlabbatch{1}.spm.spatial.smooth.data = {[tmpd,guid,'.nii,1']};
                matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2];
                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                matlabbatch{1}.spm.spatial.smooth.im = 0;
                matlabbatch{1}.spm.spatial.smooth.prefix = 's';
                spm_jobman('run',{matlabbatch});
                clear matlabbatch
                delete([tmpd,guid,'.nii']);
                fis{f}=[tmpd,'s',guid,'.nii'];
            end
        end
        n=ea_load_nii(fis{f});
        if dos || dok % cleanup temp dir
           delete(fis{f}); 
        end
        if ~exist('X','var')
            X=n.img(:);
        else
            X=[X,n.img(:)];
        end
    end