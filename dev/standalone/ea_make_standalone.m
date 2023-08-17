  function ea_make_standalone(outdir)
% Under development
%
% Compile Lead-DBS as a standalone executable using the MATLAB compiler
%   http://www.mathworks.com/products/compiler/
%
% This will generate a standalone program, which can be run
% outside MATLAB, and therefore does not use up a MATLAB licence.
%

%% Lead init...

close all
h = lead; % set paths
close(h);

outdir = fullfile(ea_getearoot,'..','lead_standalone');
if exist(outdir, 'dir'), rmdir(outdir, 's'); end % reset
mkdir(outdir);

cd(ea_getearoot); % use relative path so that a /Lead folder is created in mcr directory

%% Lead required files

% basic
required_files = {...
    fullfile('.','*.m'),...
    fullfile('.','common','ea_prefs_default.*'),...
    fullfile('.','templates','space','MNI152NLin2009bAsym','spacedef.mat'),...
    fullfile('.','templates','space','MNI152NLin2009bAsym','*.mz3'),...
    fullfile('.','templates','space','MNI152NLin2009bAsym','*.stl'),...
    fullfile('.','templates','electrode_contacts'),...
    fullfile('.','templates','electrode_models'),...
    fullfile('.','helpers'),...
    fullfile('.','ext_libs','*.m'),...
    fullfile('.','icons'),...
    fullfile('.','.version.txt')...
    };

% add more folders... (maybe optimize this?)
required_files = [required_files, {...
    fullfile('.','classes'),...
    fullfile('.','cluster'),...
    fullfile('.','connectomics'),...
    fullfile('.','classes'),...
    fullfile('.','dbshub'),...
    fullfile('.','dependency'),...
    fullfile('.','dev'),...
    fullfile('.','legacy'),...
    fullfile('.','ls'),...
    fullfile('.','predict'),...
    fullfile('.','support_scripts'),...
    fullfile('.','test'),...
    fullfile('.','toolbox'),...
    fullfile('.','tools'),...
    fullfile('.','vatmodel')...
    }];

% ext_libs (more specific inclusion instead of /ext_libs)

% required_files = [required_files, ...
%                   get_ext_libs(...
%                     'dragndrop',...
%                     'BRAINSTools',...
%                     'fsl',...
%                     'segment',...
%                     'ANTs',...
%                     'PaCER'...
%                   )];

required_files = [required_files, get_ext_libs];

%% SPM

spm fmri;
close all;

spm_standalone_part;

required_files{end+1} = spm('Dir');


%% Compilation

% check startup file. see http://www.mathworks.com/support/solutions/data/1-QXFMQ.html?1-QXFMQ
if exist('startup','file')
    warning('A startup.m has been detected in %s.\n(Use if ~isdeployed statment)',...
        fileparts(which('startup')));
end

% opts

aopts = {};
for i = 1:length(required_files)
    aopts = [aopts, {'-a', required_files{i}}];
end

Nopts = {...
    '-p',fullfile(matlabroot,'toolbox','signal'),...
    '-p',fullfile(matlabroot,'toolbox','stats'),...
    '-p',fullfile(matlabroot,'toolbox','images')...
    };

Ropts = {'-R','-singleCompThread'} ;
if spm_check_version('matlab','8.4') >= 0
    Ropts = [Ropts, {'-R','-softwareopengl'}];
end

% compile
%setenv('MCC_USE_DEPFUN','1'); % this might be necessary in case of a mcc error

mcc('-m', '-C', '-v',...
    '-o',ea_getstdaloneoname,...
    '-d',outdir,...
    '-N',Nopts{:},...
    Ropts{:},...
    aopts{:},...
    'lead.m');


%% functions

function required_ext_libs = get_ext_libs(varargin)
% get ext_libs file names to include in the compiled app.
% the inputs are the names of the required ext libs. if no inputs are
% provided then all ext libs are included

required_ext_libs = {};

if ~nargin % all ext libs
    libs_listing = dir(fullfile(ea_getearoot,'ext_libs','*'));
    libs_listing(~[libs_listing.isdir]) = []; % keep directories only
    libs_listing(strcmp('.',{libs_listing.name})) = []; % remove '.'
    libs_listing(strcmp('..',{libs_listing.name})) = []; % remove '..'
    libs_name = {libs_listing.name}; % all ext_libs.
else
    libs_name = varargin; % only the ones from input
end

% file extensions to leave out
cmp = strcmp(computer,{'PCWIN64','GLNXA64','MACI64','MACA64'});
exclude_extension = {'.exe', '.glnxa64', '.maci64', 'maca64', '.mexw64', '.mexa64', '.mexmaci64', 'mexmaca64'};
exclude_extension = exclude_extension([~cmp, ~cmp]);
exclude_extension = [exclude_extension, {'.c', '.cpp'}];

for i = 1:length(libs_name) % iterate over directories
    listing = dir(fullfile(ea_getearoot,'ext_libs',libs_name{i},'**','*')); % get all files
    listing([listing.isdir]) = []; % remove directories

    for j = 1:length(listing) % iterate over files in directory
        [filepath,name,ext] = fileparts(listing(j).name);
        if ~any(strcmp(ext,exclude_extension)) % dont include unnecesary files
            required_ext_libs{end+1} = fullfile(listing(j).folder, listing(j).name);
        end
    end
end

end

function [] = spm_standalone_part()
% this is an extract of spm_make_standalone.m

%==========================================================================
%-Static listing of SPM toolboxes
%==========================================================================
fid = fopen(fullfile(spm('Dir'),'config','spm_cfg_static_tools.m'),'wt');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
%-Get the list of toolbox directories
tbxdir = fullfile(spm('Dir'),'toolbox');
d = [tbxdir; cellstr(spm_select('FPList',tbxdir,'dir'))];
ft = {};
%-Look for '*_cfg_*.m' files in these directories
for i=1:numel(d)
    fi = spm_select('List',d{i},'.*_cfg_.*\.m$');
    if ~isempty(fi)
        ft = [ft(:); cellstr(fi)];
    end
end
%-Create code to insert toolbox config
if isempty(ft)
    ftstr = '';
else
    ft = spm_file(ft,'basename');
    ftstr = sprintf('%s ', ft{:});
end
fprintf(fid,'values = {%s};\n', ftstr);
fclose(fid);

%==========================================================================
%-Static listing of batch application initialisation files
%==========================================================================
cfg_util('dumpcfg');

%==========================================================================
% Duplicate Contents.m in Contents.txt for use in spm('Ver')
%==========================================================================
sts = copyfile(fullfile(spm('Dir'),'Contents.m'),...
    fullfile(spm('Dir'),'Contents.txt'));
if ~sts, warning('Copy of Contents.m failed.'); end

end


end
