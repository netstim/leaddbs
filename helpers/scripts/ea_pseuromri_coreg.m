function ea_pseuromri_coreg(inputlabelfile,templatecell,resolution)
% function assumes labelfile has already been roughly coregistered to
% template series
lead
if ~exist('resolution','var'); 
    resolution=0.25;
end
if ~exist('inputlabelfile','var');
    inputlabelfile='all_labels_coreg.nii';
end
if ~exist('templatecell','var');
    templatecell={['t1.nii'],['t2.nii'],['pd.nii'],'atlas.nii'};
end


labelfile='labels.nii';

copyfile([inputlabelfile],labelfile);

ea_conformspaceto(['cropped_templates',filesep,templatecell{1}],'labels.nii',0);
ea_reslice_nii(['labels.nii'],['labels.nii'],[resolution,resolution,resolution],1,0,0,[],[],1);


[pth,fname,ext]=fileparts(labelfile);
if 0
    copyfile(labelfile,fullfile(pth,['c',fname,ext]));
    labelfile=fullfile(pth,['c',fname,ext]);
    
    ea_crop_nii(labelfile,'w','nz',0,0);
    movefile(['w',labelfile],labelfile);
    
    ea_reslice_nii(labelfile,[labelfile],[0.2,0.2,0.2],1,0,0,[],[],1);
    nii=ea_load_nii(labelfile);
    nii.dt=[4,0];
    delete(labelfile);
    ea_write_nii(nii);
end



if 0 % crop templates to size of pseudo-mri
    for temp=1:length(templatecell)
        copyfile([ea_space,templatecell{temp}],['cropped_templates',filesep,templatecell{temp}]); % recopy clean version from space folder.
        ea_conformspaceto([labelfile],['cropped_templates',filesep,templatecell{temp}],1);
    end
end
opts=ea_getptopts(pth);

olabelf=ea_load_nii('labels.nii');
d=figure('Name','Progress');
for iter=1:20
    copyfile([labelfile],fullfile(pth,[fname,'_label_iter',num2str(iter-1),ext]));

    % build pseudoMRIs:
    labelf=ea_load_nii(labelfile);
    
    sim(iter)=corr(labelf.img(:),olabelf.img(:));
    set(0,'CurrentFigure',d);
    plot(sim,'bo-');
    drawnow
    olabelf=labelf;
    
    labelf.img(isnan(labelf.img))=0;

    for temp=1:length(templatecell)
        if ~exist(['r',templatecell{temp}],'file');
            copyfile(['cropped_templates',filesep,templatecell{temp}],['r',templatecell{temp}]);            
            ea_conformspaceto(labelfile,['r',templatecell{temp}],1);
        end
        tempf=ea_load_nii(['r',templatecell{temp}]);
        ixs=unique(labelf.img(:));
        ixs(ixs==0)=[];
        
        pseudof=labelf;
        pseudof.img(:)=0;
        
        for ix=ixs' % iterate through integer "colors" of the original labelfile
            pseudoixs=double(labelf.img(:)==ix).*tempf.img(:); % get all indices from template matching the pseudoMRI
            pseudof.img(labelf.img==ix)=mean(pseudoixs(~pseudoixs==0));
        end
        [tpth,tfn]=fileparts(templatecell{temp});
        pseudof.fname=fullfile(pth,['pseudo_',tfn,'_l.nii']);
        pseudof.dt=[4,0];
        ea_write_nii(pseudof);
        matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pth,['pseudo_',tfn,'_l.nii,1'])};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [0.5 0.5 0.5];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = ~strcmp('atlas.nii',templatecell{temp});
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',{matlabbatch}); clear matlabbatch
        labelcell{temp}=fullfile(pth,['spseudo_',tfn,'_l.nii']);
    end
   
    
    
    
    ea_ants_nonlinear_local(cellfun(@horzcat,repmat({['r']},size(templatecell),1),templatecell,'Uniformoutput',0),labelcell,'warped.nii',[1,1,1,2],length(templatecell),1),'MultiLabel');
    % dump output
    
    ea_ants_applytransforms(opts,{labelfile},{labelfile},0,labelcell{1},fullfile(pth,'warpedComposite.h5'),'MultiLabel');

end









