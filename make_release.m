function make_release(outdir)
% make release

if ~strcmp(outdir(end),filesep)
    outdir=[outdir,filesep];
end

mkdir([outdir,'lead']);
copyfile([fileparts(which('lead'))],[outdir,'lead']);

% reset prefs:
movefile([outdir,'lead',filesep,'ea_prefs_public.m'],[outdir,'lead',filesep,'ea_prefs.m']);
delete([outdir,'lead',filesep,'ea_prefs.mat']);

% delete git:
%rmdir([outdir,'lead',filesep,'.git'],'s');

% delete atlases:

rmdir([outdir,'lead',filesep,'atlases'],'s');
mkdir([outdir,'lead',filesep,'atlases']);
mkdir([outdir,'lead',filesep,'atlases',filesep,'Put atlases here']);
mkdir([outdir,'lead',filesep,'atlases',filesep,'atlasset_1']);
mkdir([outdir,'lead',filesep,'atlases',filesep,'atlasset_1',filesep,'lh']);
mkdir([outdir,'lead',filesep,'atlases',filesep,'atlasset_1',filesep,'rh']);
mkdir([outdir,'lead',filesep,'atlases',filesep,'atlasset_1',filesep,'mixed']);
mkdir([outdir,'lead',filesep,'atlases',filesep,'atlasset_1',filesep,'midline']);


% delete cfg:

rmdir([outdir,'lead',filesep,'cfg'],'s');

% delete trajvectors:

delete([outdir,'lead',filesep,'trajvectors.mat']);

% delete make_release:

delete([outdir,'lead',filesep,'make_release.m']);


% delete DARTEL-Part:

rmdir([outdir,'lead',filesep,'templates',filesep,'dartel'],'s');
delete([outdir,'lead',filesep,'ea_normalize_spmda*']);

% remove Gibbshighest:
delete([outdir,'lead',filesep,'fibers',filesep,'gibbsconnectome.mat']);
