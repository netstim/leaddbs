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
    outfn=[ea_space,'atlas'];
end

[pth,fn,ext]=fileparts(outfn);
outfn=fullfile(pth,ea_stripext(fn));

[infn,ext]=ea_niigz([ea_space([],'atlases'),atlassetname,filesep,'gm_mask.nii']);
copyfile(infn,[outfn,ext]);
switch ext
    case '.nii.gz'
        gunzip([outfn,ext]);
        ea_delete([outfn,ext]);
end

ea_conformspaceto([ea_space,spacedef.templates{1},'.nii'],[outfn,'.nii'],1);
