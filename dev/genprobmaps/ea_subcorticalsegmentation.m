function ea_subcorticalsegmentation(options)

directory=[options.root,options.patientname,filesep];

if exist([directory,'mmprocessed.mat'],'file');
    return
end

%% find multi files

if ~exist([directory,'mmpreprocessed.mat'],'file');
    ea_preparesubcorticalsegmentation(options);
end


mults=dir([directory,'sbbmult*.nii']);
if isempty(mults)
    ea_error('Please place additional MR modality acquisitions into patient folder (called mult1.nii, mult2.nii and so on.');
end
srcs{1}=[directory,'sbb',options.prefs.prenii_unnormalized];
for m=1:length(mults)
    srcs{1+m}=[directory,mults(m).name];
end



%% start multimodal subcortical segmentation algorithm.
structures={'STN','Pallidum','Ruber'};
%structures={'wm','gm','csf'};

genmaps=1;
wtamaps=1;
cleanup=0;

if genmaps
    for s=1:length(structures)
        pts=ea_readcsv([structures{s},'.fcsv']);
        
        [~,tempfile]=ea_whichnormmethod(directory);
        V=spm_vol(tempfile);
        pts=V.mat\pts;
        pts=ea_map_coords(pts,tempfile,[directory,'y_ea_normparams.nii'],[directory,options.prefs.prenii_unnormalized]);
        label=structures{s};
        
        ea_generate_probmaps(pts(1:3,:)',label,srcs,directory);
    end
end

if wtamaps
    
    dfactor=1.2;
    
    for s=1:length(structures)
        S{s}=ea_load_nii([directory,'s',structures{s},'_secondlevel.nii']);
        A(:,:,:,s)=S{s}.img;
        S{s}.img(:)=nan;
    end
    
    ea_dispercent(0,'WTAing');
    for xx=1:size(S{s}.img,1)
        for yy=1:size(S{s}.img,2)
            for zz=1:size(S{s}.img,3)
                vals=squeeze(A(xx,yy,zz,:));
                if any(vals>0.3)
                    [v,ix]=sort(vals);
                    if v(end)>dfactor*v(end-1)
                        S{ix(end)}.img(xx,yy,zz)=1;
                    end
                end
            end
        end
        ea_dispercent(xx/size(S{s}.img,1));
    end
    ea_dispercent(1,'end');
    % put in atlases dir:
    mkdir([directory,'atlases']);
    mkdir([directory,'atlases',filesep,'native']);
    atldir=[directory,'atlases',filesep,'native',filesep,options.patientname];
    mkdir(atldir);
    mkdir([atldir,filesep,'lh']);
    mkdir([atldir,filesep,'rh']);
    for s=1:length(structures)
        S{s}.fname=[directory,structures{s},'_wta.nii'];
        spm_write_vol(S{s},S{s}.img);
        ea_largestcomponent_nii([directory,structures{s},'_wta.nii'],'hemispheres');
        
        movefile([directory,structures{s},'_wta_lh.nii'],[atldir,filesep,'lh',filesep,structures{s},'.nii']);
        movefile([directory,structures{s},'_wta_rh.nii'],[atldir,filesep,'rh',filesep,structures{s},'.nii']);
        ea_crop_nii([atldir,filesep,'lh',filesep,structures{s},'.nii']);
        ea_crop_nii([atldir,filesep,'rh',filesep,structures{s},'.nii']);
    end
end


if cleanup
    for s=1:length(structures)
        try       delete([directory,'',structures{s},'_secondlevel.nii']); end
        try       delete([directory,'',structures{s},'_firstlevel.nii']); end
        try       delete([directory,'',structures{s},'_wta.nii']); end
        try       delete([directory,'s',structures{s},'_secondlevel.nii']); end
    end
end

mmprocessed=1;
save([directory,'mmprocessed.mat'],'mmprocessed');


function ea_preparesubcorticalsegmentation(options)
directory=[options.root,options.patientname,filesep];

mults=dir([directory,'mult*.nii']);
if isempty(mults)
    ea_error('Please place additional MR modality acquisitions into patient folder (called multi1.nii, multi2.nii and so on.');
end
ea_normalize_reslicepretra(options);

[whichnormmethod]=ea_whichnormmethod(directory);
% warp bb into folder.

switch whichnormmethod
    case ea_getantsnormfuns
        ea_ants_applytransforms(options,{[options.earoot,'templates',filesep,'bb.nii']},{[directory,'wbb.nii']},1);
    otherwise
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_normparams.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.push.fnames = {[options.earoot,'templates',filesep,'bb.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
        matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {directory};
        matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {[directory,options.prenii_unnormalized]};
        matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
        matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
        cfg_util('run',{matlabbatch});
        clear matlabbatch
end

ea_reslice_nii([directory,'wbb.nii'],[directory,'wbb.nii'],[1,1,1]);

ea_crop_nii([directory,'wbb.nii']);
ea_reslice_nii([directory,'wbb.nii'],[directory,'wbb.nii'],[0.22,0.22,0.22]);



% anat

matlabbatch{1}.spm.util.imcalc.input = {[directory,'wbb.nii']
    [directory,options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.util.imcalc.output = ['bb',options.prefs.prenii_unnormalized];
matlabbatch{1}.spm.util.imcalc.outdir = {directory};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
cfg_util('run',{matlabbatch});
clear matlabbatch

matlabbatch{1}.spm.spatial.smooth.data = {[directory,'bb',options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
cfg_util('run',{matlabbatch});
clear matlabbatch

delete([directory,'bb',options.prefs.prenii_unnormalized]);


% mults

for m=1:length(mults)
    try movefile([directory,'raw_',mults(m).name],[directory,mults(m).name]); end
    copyfile([directory,mults(m).name],[directory,'raw_',mults(m).name]);
    switch options.coregmr.method
        case 1 % SPM
            ea_docoreg_spm([directory,mults(m).name],[directory,options.prefs.prenii_unnormalized],'nmi',1)
        case 2 % ANTs
            ea_ants([directory,options.prefs.prenii_unnormalized],...
                [directory,mults(m).name],...
                [directory,mults(m).name],0);
        case 3 % BRAINSFit
            ea_brainsfit([directory,options.prefs.prenii_unnormalized],...
                [directory,mults(m).name],...
                [directory,mults(m).name],0);
        case 4 % Hybrid SPM -> ANTs
            ea_docoreg_spm([directory,mults(m).name],[directory,options.prefs.prenii_unnormalized],'nmi',0)
            ea_ants([directory,options.prefs.prenii_unnormalized],...
                [directory,mults(m).name],...
                [directory,mults(m).name],0);
        case 5 % Hybrid SPM -> Brainsfit
            ea_docoreg_spm([directory,mults(m).name],[directory,options.prefs.prenii_unnormalized],'nmi',0)
            ea_brainsfit([directory,options.prefs.prenii_unnormalized],...
                [directory,mults(m).name],...
                [directory,mults(m).name],0);
    end
    
    
    matlabbatch{1}.spm.util.imcalc.input = {[directory,'wbb.nii']
        [directory,mults(m).name]};
    matlabbatch{1}.spm.util.imcalc.output = ['bb',mults(m).name];
    matlabbatch{1}.spm.util.imcalc.outdir = {directory};
    matlabbatch{1}.spm.util.imcalc.expression = 'i2';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    cfg_util('run',{matlabbatch});
    clear matlabbatch
    
    matlabbatch{1}.spm.spatial.smooth.data = {[directory,'bb',mults(m).name]};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    cfg_util('run',{matlabbatch});
    clear matlabbatch
    
    delete([directory,'bb',mults(m).name]);
end

delete([directory,'wbb.nii']);
mmpreprocessed=1;
save([directory,'mmpreprocessed.mat'],'mmpreprocessed');


function c=ea_readcsv(pth)
fid=fopen(pth);
C=textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %s %s %s','commentStyle', '#','delimiter', ',');
fclose(fid);
for coord=1:length(C{1})
    c(:,coord)=[C{2}(coord);C{3}(coord);C{4}(coord);1];
end