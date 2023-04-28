function ea_genauxspace(~,~,handles)
options.prefs=ea_prefs('');
if exist([ea_space,'spacedef.mat'],'file')
    load([ea_space,'spacedef.mat'])
else
    spacedef=ea_gendefspacedef; % generate default spacedefiniton
end

%% 3. Create atlas.nii from selected atlas:
if ischar(handles)
    atlassetname=handles;
else
    atlassetname=get(handles.atlassetpopup,'String');
    if ~isempty(atlassetname)
        atlassetname=atlassetname{get(handles.atlassetpopup,'Value')};  
    end
end

ea_flattenatlas(atlassetname);

%% 4. Create TPM, DARTEL and SHOOT templates:
ea_create_tpm_darteltemplate('mute')

%% 5. binarize c1 and c2 masks that come out of this as a by-product:
ea_binmasks

%% 6. Create Wires:
if ~exist([ea_space,'wires.mat'],'file')
    ea_genwires
end

%% Warp remaining assets from MNI 2009b NLIN

if ~exist([ea_space,'subcortical'], 'dir')
    mkdir([ea_space,'subcortical']);
end
ea_importspaceassets([],[],'MNI152NLin2009bAsym','custom',[ea_getearoot,'templates',filesep,'space',filesep,'MNI152NLin2009bAsym',filesep,'subcortical',filesep,'secondstepmask.nii'],[ea_space(options,'subcortical'),'secondstepmask.nii'])
ea_importspaceassets([],[],'MNI152NLin2009bAsym','custom',[ea_getearoot,'templates',filesep,'space',filesep,'MNI152NLin2009bAsym',filesep,'subcortical',filesep,'thirdstepmask.nii'],[ea_space(options,'subcortical'),'thirdstepmask.nii'])


ea_importspaceassets([],[],'MNI152NLin2009bAsym','custom',[ea_getearoot,'templates',filesep,'space',filesep,'MNI152NLin2009bAsym',filesep,'bb.nii'],[ea_space,'bb.nii'])
ea_importspaceassets([],[],'MNI152NLin2009bAsym','custom',[ea_getearoot,'templates',filesep,'space',filesep,'MNI152NLin2009bAsym',filesep,'t1.nii'],[ea_space,'tempt.nii'])
%ea_importspaceassets([],[],'MNI152NLin2009bAsym','delete') % remove warp directory again.

nii=ea_load_nii([ea_space,'bb.nii']);
nii.img=nii.img>0.1;
ea_write_nii(nii);
ea_crop_nii([ea_space,'bb.nii']);


function ea_binmasks
nii=ea_load_nii([ea_space,'c1mask.nii']);
delete([ea_space,'c1mask.nii']);
nii.img=nii.img>0.5;
nii.dt(1) = 2;
ea_write_nii(nii);
nii=ea_load_nii([ea_space,'c2mask.nii']);
delete([ea_space,'c2mask.nii']);
nii.img=nii.img>0.5;
nii.dt(1) = 2;
ea_write_nii(nii);
