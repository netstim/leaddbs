function make_release(varargin)
% make release

if ~nargin % default output.
    
    outdir='/PA/Neuro/_projects/lead/release';
    disp(['Outputting to /PA/Neuro/_projects/lead/release']);
else
    outdir=varargin{1};
end

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

leave_atlases={'ATAG_Linear','ATAG_Nonlinear'};

atls=dir([outdir,'lead',filesep,'atlases']);

for atl=1:length(atls);
    if ~ismember(atls(atl).name,leave_atlases) && ~strcmp(atls(atl).name,'.') && ~strcmp(atls(atl).name,'..') && ~strcmp(atls(atl).name,'.DS_Store')
try        
        rmdir([outdir,'lead',filesep,'atlases',filesep,atls(atl).name],'s');
catch
    disp(['Couldnt delete atlases',filesep,atls(atl).name,'.']);
end
    end
end
%mkdir([outdir,'lead',filesep,'atlases']);


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
