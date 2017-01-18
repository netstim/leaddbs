function ea_genauxspace(~,~,handles)


options.prefs=ea_prefs('');
if exist([ea_space,'ea_space_def.mat'],'file')
    load([ea_space,'ea_space_def.mat'])
else
    spacedef=ea_gendefspacedef; % generate default spacedefiniton
end
%% 1. generate bb.nii

copyfile([ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'bb.nii'],[ea_space,'bb.nii']);
copyfile([ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'t1.nii'],[ea_space,'tempt.nii']);
copyfile([ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'subcortical',filesep,'secondstepmask.nii'],[ea_space,'secondstepmask.nii']);
copyfile([ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'subcortical',filesep,'thirdstepmask.nii'],[ea_space,'thirdstepmask.nii']);


matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ea_space,spacedef.templates{1},'.nii,1']}; % for now assume there is a T1 at least in that new space..
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[ea_space,'tempt.nii,1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[ea_space,'bb.nii,1']
    [ea_space,'secondstepmask.nii']
    [ea_space,'thirdstepmask.nii']};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',{matlabbatch});
clear matlabbatch
delete([ea_space,'rbb.nii']);
delete([ea_space,'rtempt.nii']);
delete([ea_space,'tempt.nii']);


matlabbatch{1}.spm.util.imcalc.input = {
                                        [ea_space,'bb.nii,1']
                                        [ea_space,spacedef.templates{1},'.nii,1']
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'bb.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {ea_space};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',{matlabbatch});
clear matlabbatch

%% 2. Create "subcortical" folder:
mkdir([ea_space,'subcortical']);
movefile([ea_space,'rsecondstepmask.nii'],[ea_space,'subcortical',filesep,'secondstepmask.nii']);
movefile([ea_space,'rthirdstepmask.nii'],[ea_space,'subcortical',filesep,'thirdstepmask.nii']);
delete([ea_space,'thirdstepmask.nii']); delete([ea_space,'secondstepmask.nii']);

%% 3. Create atlas.nii from selected atlas:
if ischar(handles)
    atlassetname=handles;
else
    atlassetname=get(handles.atlassetpopup,'String');
    if ~isempty(atlassetname)
        atlassetname=atlassetname{get(handles.atlassetpopup,'Value')};
        ea_flattenatlas(atlassetname);
    end
end

%% 4. Create TPM, DARTEL and SHOOT templates:
ea_create_tpm_darteltemplate('mute')

%% 5. binarize c1 and c2 masks that come out of this as a by-product:
ea_binmasks

%% 6. Create Wires:
ea_genwires



function ea_binmasks
nii=ea_load_nii([ea_space,'c1mask.nii']);
delete([ea_space,'c1mask.nii']);
nii.img=nii.img>0.5;
nii.dt=[2,0];
ea_write_nii(nii);
nii=ea_load_nii([ea_space,'c2mask.nii']);
delete([ea_space,'c2mask.nii']);
nii.img=nii.img>0.5;
nii.dt=[2,0];
ea_write_nii(nii);
