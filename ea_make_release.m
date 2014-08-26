function ea_make_release(varargin)
% make release - give increase of code and content version number.
 
outdir='/PA/Neuro/_projects/lead/release';
disp(['Outputting to /PA/Neuro/_projects/lead/release']);
    
if ~nargin % default output.
    inc_code=0.01;
    inc_cont=0;
elseif nargin==1
    inc_code=varargin{1};
    inc_cont=0;
elseif nargin==2
    inc_code=varargin{1};
    inc_cont=varargin{2};
end

if ~strcmp(outdir(end),filesep)
    outdir=[outdir,filesep];
end

try rmdir([outdir,'lead_dbs'],'s'); end
mkdir([outdir,'lead_dbs']);
copyfile([fileparts(which('lead'))],[outdir,'lead_dbs']);

% reset prefs:
movefile([outdir,'lead_dbs',filesep,'ea_prefs_public.m'],[outdir,'lead_dbs',filesep,'ea_prefs.m']);
delete([outdir,'lead_dbs',filesep,'ea_prefs.mat']);
delete([outdir,'lead_dbs',filesep,'ea_ui.mat']);


% delete atlases:

leave_atlases={'ATAG_Linear (Keuken 2014)','ATAG_Nonlinear (Keuken 2014)','ATAG_STN (Forstmann 2012 & Keuken 2013)'};

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


% delete ui:

delete([outdir,'lead_dbs',filesep,'ea_ui.mat']);


% delete make_release:

delete([outdir,'lead_dbs',filesep,'make_release.m']);


% delete dev folder:

rmdir([outdir,'lead_dbs',filesep,'dev'],'s');


% delete DARTEL-Templates (can be generated the first time they are used):

delete([outdir,'lead_dbs',filesep,'templates',filesep,'dartel',filesep,'dartelmni_*.nii']);
%delete([outdir,'lead_dbs',filesep,'ea_normalize_spmda*']);



%% zip all_content:
addfolders={[outdir,'lead_dbs']};
zip([outdir,'lead_full_release.zip'],addfolders);





%% upload to FTP:
disp('Connecting to FTP-Server...');
mw = ftp('www.andreas-horn.de','www.andreas-horn.de','andiANDI$1');
disp('Changing Dir.');
cd(mw,'leaddbs/release');
disp('Uploading full release.');
mput(mw, [outdir,'lead_full_release.zip']);
disp('Done.');
close(mw);


%% zip addcontent:
addfolders={[outdir,'lead_dbs',filesep,'templates'],[outdir,'lead_dbs',filesep,'fibers'],[outdir,'lead_dbs',filesep,'atlases']};
zip([outdir,'lead_content.zip'],addfolders);





%% upload to FTP:
disp('Connecting to FTP-Server...');
mw = ftp('www.andreas-horn.de','www.andreas-horn.de','andiANDI$1');
disp('Changing Dir.');
cd(mw,'leaddbs/release');
disp('Uploading content.');
mput(mw, [outdir,'lead_content.zip']);
disp('Done.');
close(mw);


%% remove addcontent:
for fi=1:length(addfolders)
    rmdir(addfolders{fi},'s');
end


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


%% update version:

v=ea_getvsn('local');
v=v+[inc_code;inc_cont];
delete([fileparts(which('lead')),filesep,'.version.txt']);
fileID = fopen([fileparts(which('lead')),filesep,'.version.txt'],'w');
fprintf(fileID,'%6.3f\n',v);
fclose(fileID);


%% upload version to FTP:
disp('Connecting to FTP-Server...');
mw = ftp('www.andreas-horn.de','www.andreas-horn.de','andiANDI$1');
disp('Changing Dir.');
cd(mw,'leaddbs/release');
disp('Uploading Version file.');
mput(mw, [fileparts(which('lead')),filesep,'.version.txt']);
disp('Done.');
close(mw);


%% cleanup:
%delete([outdir,'lead_content.zip']);
%delete([outdir,'lead_dbs.zip']);
%delete([outdir,'lead_full_release.zip']);

