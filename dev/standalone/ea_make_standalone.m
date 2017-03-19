function ea_make_standalone(outdir)
% Compile Lead-DBS as a standalone executable using the MATLAB compiler
%   http://www.mathworks.com/products/compiler/
%
% This will generate a standalone program, which can be run
% outside MATLAB, and therefore does not use up a MATLAB licence.
%
%% Lead part...
close all
lead % set paths
close all
outdir = fullfile(ea_getearoot,'../lead_standalone');
mkdir(outdir);
%% make a copy of required lead DBS folders:
tof=[ea_getearoot,'standalone_export',filesep,ea_getstdaloneoname,filesep];

if exist(fileparts(fileparts(tof)),'dir')
    rmdir(fileparts(fileparts(tof)),'s');
end
mkdir(tof);
from=ea_getearoot;
copyfile([from,'*.m'],tof);
copyfile([from,'*.fig'],tof);
movefile([tof,'lead.m'],[fileparts(fileparts(tof)),filesep,'lead.m']);
movefile([tof,'lead.fig'],[fileparts(fileparts(tof)),filesep,'lead.fig']);

% helpers
copyfile([from,'helpers'],[tof,'helpers']);
% ext. libs
copyfile([from,'ext_libs'],[tof,'ext_libs']);
delete([tof,'ext_libs',filesep,'*.c']);
delete([tof,'ext_libs',filesep,'*.cpp']);
delete([tof,'ext_libs',filesep,'*/*.c']);
delete([tof,'ext_libs',filesep,'*/*.cpp']);
delete([tof,'ext_libs',filesep,'*/*/*.c']);
delete([tof,'ext_libs',filesep,'*/*/*.cpp']);
delete([tof,'ext_libs',filesep,'*/*/*/*.c']);
delete([tof,'ext_libs',filesep,'*/*/*/*.cpp']);
% connectomics
copyfile([from,'connectomics'],[tof,'connectomics']);

% icons
copyfile([from,'icons'],[tof,'icons']);

% vatmodel
copyfile([from,'vatmodel'],[tof,'vatmodel']);

% tools
copyfile([from,'tools'],[tof,'tools']);

% templates
mkdir([tof,'templates',filesep,'space',filesep]);
copyfile([from,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM'],[tof,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM']);
rmdir([tof,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'atlases'],'s');
mkdir([tof,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'atlases']);
copyfile([from,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'atlases',filesep,'DISTAL (Ewert 2016)'],[tof,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'atlases',filesep,'DISTAL (Ewert 2016)']);
copyfile([from,'templates',filesep,'electrode_contacts'],[tof,'templates',filesep,'electrode_contacts']);
copyfile([from,'templates',filesep,'electrode_models'],[tof,'templates',filesep,'electrode_models']);

delete([tof,'common',filesep,'ea_recentpatients.mat']);

rmpath(genpath(from));
addpath(genpath(tof));
cd(tof);

% for now delete some ext_libs that are not needed:
rmdir([tof,'ext_libs',filesep,'surfice'],'s');
rmdir([tof,'ext_libs',filesep,'dsi_studio'],'s');




%% SPM part...
spm fmri
close all



%__________________________________________________________________________

%--------------------------------------------------------------------------
% see http://www.mathworks.com/support/solutions/data/1-QXFMQ.html?1-QXFMQ
if exist('startup','file')
    warning('A startup.m has been detected in %s.\n',...
        fileparts(which('startup')));
end
try rmdir(fullfile(ea_getearoot,'../../lead_standalone'),'s'); end


%==========================================================================
%-Static listing of SPM toolboxes
%==========================================================================
fid = fopen(fullfile(spm('dir'),'config','spm_cfg_static_tools.m'),'wt');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
% create code to insert toolbox config
%-Toolbox autodetection
%-Get the list of toolbox directories
tbxdir = fullfile(spm('Dir'),'toolbox');
d  = dir(tbxdir); d = {d([d.isdir]).name};
dd = regexp(d,'^\.');
%(Beware, regexp returns an array if input cell array is of dim 0 or 1)
if ~iscell(dd), dd = {dd}; end
d  = {'' d{cellfun('isempty',dd)}};
ft = {};
%-Look for '*_cfg_*.m' files in these directories
for i=1:length(d)
    d2 = fullfile(tbxdir,d{i});
    di = dir(d2); di = {di(~[di.isdir]).name};
    f2 = regexp(di,'.*_cfg_.*\.m$');
    if ~iscell(f2), f2 = {f2}; end
    fi = di(~cellfun('isempty',f2));
    if ~isempty(fi)
        ft = [ft(:); fi(:)];
    end
end
if ~isempty(ft)
    if isempty(ft)
        ftstr = '';
    else
        ft = cellfun(@(cft)strtok(cft,'.'),ft,'UniformOutput',false);
        ftstr  = sprintf('%s ', ft{:});
    end
    fprintf(fid,'values = {%s};\n', ftstr);
end
fclose(fid);

%==========================================================================
%-Static listing of batch application initialisation files
%==========================================================================
cfg_util('dumpcfg');

%==========================================================================
%-Duplicate Contents.m in Contents.txt for use in spm('Ver')
%==========================================================================
sts = copyfile(fullfile(spm('Dir'),'Contents.m'),...
               fullfile(spm('Dir'),'Contents.txt'));
if ~sts, warning('Copy of Contents.m failed.'); end

%==========================================================================
%-Compilation
%==========================================================================
opts = {'-p',fullfile(matlabroot,'toolbox','signal')};
mcc('-m', '-C', '-v',...
    '-o',ea_getstdaloneoname,...
    '-d',outdir,...
    '-N',opts{:},...
    '-R','-singleCompThread',...
    '-a',fileparts(ea_getearoot),...
    'lead.m');
rmpath(genpath(tof));
rmdir(tof,'s');
addpath(from);
