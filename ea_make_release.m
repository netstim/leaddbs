function ea_make_release(varargin)
% make release - give increase of code and content version number.
% ea_make_release(outdir,inc_code);

outdir='/PA/Neuro/_projects/lead/release';
    
if ~nargin % default output.
    inc_code=0.001;
elseif nargin==1
        inc_code=0.001;
    if ~isempty(varargin{1})
    outdir=varargin{1};
    end
elseif nargin==2
if ~isempty(varargin{1})
    outdir=varargin{1};
end
inc_code=varargin{2};
elseif nargin==3
if ~isempty(varargin{1})
    outdir=varargin{1};
end
inc_code=varargin{2};
end

disp(['Outputting to ',outdir,'.']);


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

leave_atlases={'CFA Subcortical Shape Atlas (Qiu 2010)','ATAG_Linear (Keuken 2014)','ATAG_Nonlinear (Keuken 2014)','ATAG_STN (Forstmann 2012 & Keuken 2013)','STN-Subdivisions (Accolla 2014)','BGHAT (Prodoehl 2008)'};

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

% delete labelings:

leave_atlases={'aal.nii','aal.txt'};

atls=dir([outdir,'lead_dbs',filesep,'templates',filesep,'labeling']);

for atl=1:length(atls);
    if ~ismember(atls(atl).name,leave_atlases) && ~strcmp(atls(atl).name,'.') && ~strcmp(atls(atl).name,'..') && ~strcmp(atls(atl).name,'.DS_Store')
try        
        rm([outdir,'lead_dbs',filesep,'atlases',filesep,atls(atl).name]);
catch
    disp(['Couldnt delete atlases',filesep,atls(atl).name,'.']);
end
    end
end


%mkdir([outdir,'lead_dbs',filesep,'atlases']);


% delete cfg:

rmdir([outdir,'lead_dbs',filesep,'cfg'],'s');


% delete Manual:

rmdir([outdir,'lead_dbs',filesep,'Lead_manual'],'s');


% delete trajvectors:

delete([outdir,'lead_dbs',filesep,'trajvectors.mat']);


% delete modeldti:
delete([outdir,'lead_dbs',filesep,'ea_genvat_modeldti.m']);

% delete ui:

delete([outdir,'lead_dbs',filesep,'ea_ui.mat']);

% delete TPM.nii:

delete([outdir,'lead_dbs',filesep,'templates',filesep,'TPM.nii']);


% delete make_release:

delete([outdir,'lead_dbs',filesep,'ea_make_release.m']);


% delete dev folder:

rmdir([outdir,'lead_dbs',filesep,'dev'],'s');

% delete fibers:

delete([outdir,'lead_dbs',filesep,'fibers',filesep,'*']);


% delete tmp folder:

try rmdir([outdir,'lead_dbs',filesep,'tmp'],'s'); end

% delete DARTEL-Templates (can be generated the first time they are used):

delete([outdir,'lead_dbs',filesep,'templates',filesep,'dartel',filesep,'dartelmni_*.nii']);
%delete([outdir,'lead_dbs',filesep,'ea_normalize_spmda*']);


%% update version of release:

v=ea_getvsn('local');
v=v+inc_code;
%delete([outdir,'lead_dbs',filesep,'.version.txt']);
fileID = fopen([outdir,'lead_dbs',filesep,'.version.txt'],'w');
fprintf(fileID,'%6.3f\n',v);
fclose(fileID);




%% zip all_content:
addfolders={[outdir,'lead_dbs']};
zip([outdir,'lead_dbs.zip'],addfolders);





%% upload to FTP:
disp('Connecting to FTP-Server...');
mw = ftp('www.andreas-horn.de','www.andreas-horn.de','andiANDI$1');
disp('Changing Dir.');
cd(mw,'leaddbs/release');
disp('Uploading full release.');
mput(mw, [outdir,'lead_dbs.zip']);
disp('Done.');
close(mw);

delete([outdir,'lead_dbs.zip']);



%% update version (locally):

v=ea_getvsn('local');
v=v+[inc_code];
delete([fileparts(which('lead')),filesep,'.version.txt']);
fileID = fopen([fileparts(which('lead')),filesep,'.version.txt'],'w');
fprintf(fileID,'%6.3f\n',v);
fclose(fileID);

%% upload version textfile to FTP:
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

