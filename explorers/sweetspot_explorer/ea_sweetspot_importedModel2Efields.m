function [S]=ea_sweetspot_importedModel2Efields(obj, filepath)                
    %load the file
    S = load(filepath);
    disp('Need to put imported model to the same space as the Efields.');
    outdir = fullfile(fileparts(obj.leadgroup), 'sweetspots', obj.ID, filesep);
    ea_mkdir(outdir);
    
    sidesuffices = {'_r', '_l'};

    try 
        hemilateral=evalin('base','hemilateral');
    end 

    if exist('hemilateral','var')
        space=obj.results.space{1,1};
        space.fname=[outdir, 'efield_space.nii'];
        ea_write_nii(space)

        toreslice=S.space{1,hemilateral};
        toreslice.img(:)=S.model_vals{1,hemilateral};
        toreslice.fname=[outdir, 'model_space', sidesuffices{hemilateral}, '.nii'];
        ea_write_nii(toreslice)

        ea_conformspaceto([outdir,'efield_space.nii'], ...
            [outdir, 'model_space', sidesuffices{hemilateral}, '.nii'], 0);
        nii=ea_load_nii([outdir, 'model_space', sidesuffices{hemilateral}, '.nii']);
        % put 0 to NaN
        nii.img(nii.img==0)=NaN;
%         ix=nii.img==0;
%         nii.img(ix)=NaN;
        ea_write_nii(nii)       

        vals{hemilateral}=nii.img(:);
    
    else

        for side=1:numel(obj.results.space)
            space=obj.results.space{1,side};
            space.fname=[outdir, 'efield_space', sidesuffices{side}, '.nii'];
            ea_write_nii(space)
    
            toreslice=S.space{1,side};
            toreslice.img(:)=S.model_vals{1,side};
            toreslice.fname=[outdir, 'model_space', sidesuffices{side}, '.nii'];
            ea_write_nii(toreslice)
    
            ea_conformspaceto([outdir,'efield_space',sidesuffices{side},'.nii'], ...
                [outdir, 'model_space', sidesuffices{side}, '.nii'], 0);
            nii=ea_load_nii([outdir, 'model_space', sidesuffices{side}, '.nii']);
            nii.img(nii.img==0)=NaN;
            ea_write_nii(nii)
            vals{side}=nii.img(:);
        end
    end
    S.space=space;
    S.model_vals=vals;

%     ea_delete(outdir);

end
