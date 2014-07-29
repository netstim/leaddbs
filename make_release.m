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

rmdir([outdir,'lead_dbs'],'s');
mkdir([outdir,'lead_dbs']);
copyfile([fileparts(which('lead'))],[outdir,'lead_dbs']);

% reset prefs:
movefile([outdir,'lead_dbs',filesep,'ea_prefs_public.m'],[outdir,'lead_dbs',filesep,'ea_prefs.m']);
delete([outdir,'lead_dbs',filesep,'ea_prefs.mat']);

% delete git:
%rmdir([outdir,'lead_dbs',filesep,'.git'],'s');

% delete TPM file
delete([outdir,'lead_dbs',filesep,'templates',filesep,'TPM.nii']);

% delete atlases:

leave_atlases={'ATAG_Linear (Keuken 2014)','ATAG_Nonlinear (Keuken 2014)'};

atls=dir([outdir,'lead_dbs',filesep,'atlases']);

for atl=1:length(atls);
    if ~ismember(atls(atl).name,leave_atlases) && ~strcmp(atls(atl).name,'.') && ~strcmp(atls(atl).name,'..') && ~strcmp(atls(atl).name,'.DS_Store')
try        
        rmdir([outdir,'lead_dbs',filesep,'atlases',filesep,atls(atl).name],'s');
catch
    disp(['Couldnt delete atlases',filesep,atls(atl).name,'.']);
end
    end
end
%mkdir([outdir,'lead_dbs',filesep,'atlases']);


% delete cfg:

rmdir([outdir,'lead_dbs',filesep,'cfg'],'s');

% delete trajvectors:

delete([outdir,'lead_dbs',filesep,'trajvectors.mat']);

% delete make_release:

delete([outdir,'lead_dbs',filesep,'make_release.m']);


% delete dev folder:

rmdir([outdir,'lead_dbs',filesep,'dev'],'s');


% delete DARTEL-Templates (can be generated the first time they are used):

delete([outdir,'lead_dbs',filesep,'templates',filesep,'dartel',filesep,'dartelmni_*.nii']);
%delete([outdir,'lead_dbs',filesep,'ea_normalize_spmda*']);

% remove Gibbshighest:
delete([outdir,'lead_dbs',filesep,'fibers',filesep,'gibbsconnectome.mat']);



%% zip release:
zip([outdir,'lead_dbs.zip'],[outdir,'lead_dbs']);
rmdir([outdir,'lead_dbs'],'s');

%% upload to FTP:
disp('Connecting to FTP-Server...');
mw = ftp('www.andreas-horn.de','www.andreas-horn.de','andiANDI$1');
disp('Changing Dir.');
cd(mw,'leaddbs/release');
disp('Uploading release.');
mput(mw, [outdir,'lead_dbs.zip']);
disp('Done.');
close(mw);
