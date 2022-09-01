function ea_flattenatlas(atlassetname,outfn)

if isempty(atlassetname)
    return
end
load([ea_space,'spacedef.mat']);
if exist([ea_space([],'atlases'),atlassetname,filesep,'atlas_index.mat'],'file')
    load([ea_space([],'atlases'),atlassetname,filesep,'atlas_index.mat']);
else
    
    options.atlasset=atlassetname;
    options.atl.can=1;
    options.atl.ptnative=0;
    options.prefs=ea_prefs('');
    ea_genatlastable([],ea_space([],'atlases'),options)
    load([ea_space([],'atlases'),atlassetname,filesep,'atlas_index.mat']);

end

if ~exist('outfn','var')
    outfn=[ea_space,'atlas.nii'];
end
copyfile([ea_space([],'atlases'),atlassetname,filesep,'gm_mask.nii'],outfn);

ea_conformspaceto([ea_space,spacedef.templates{1},'.nii'],outfn,1);
