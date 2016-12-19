function varargout = migrate_fibertools_jobs(varargin)

%MIGRATE_FIBERTOOLS_JOBS Utility to migrate batch jobs from vY to vX
%   In release vX of DTI&Fibertools, some details of the batch
%   configuration have been changed compared to vY. While these changes
%   should make it easier to create batch jobs that are both readable and
%   concise, they introduce incompatible changes to existing jobs. Loading 
%   vY jobs directly in the vX batch configuration will lead to information
%   loss.
%   MIGRATE_FIBERTOOLS_JOBS must be used to adapt existing jobs to the new
%   configuration before these jobs can be loaded and run. Wherever
%   possible, existing information will be transferred into the new job 
%   structure. The tool will go through the batch structures supplied and
%   look for any entries that refer to vY modules of DTI&Fibertools. These
%   will be replaced with vX compatible syntax. If called without an input
%   argument, the tool will ask for a list of job files to be migrated.
%   The migrated jobs will be saved in the same place and format as the
%   original jobs, but with a '_vX' appended to the original filename.
%   It is safe to run the migration multiple times and also on jobs which
%   do not have this compatibility problem.
%   [vXjob1 ...] = MIGRATE_FIBERTOOLS_JOBS(vYjob1, ...) takes a list of
%   jobs (this can be a mixture of job filenames and job variables) and
%   runs the migration on all of them. Jobs which are read from files will
%   be saved as files as described above. Jobs which are passed as job
%   variables will be returned as variables.

% (C) 

%#ok<*TRYNC>

%% Collect jobs
if nargin == 0
    vYjobs = cfg_getfile(Inf, 'batch', 'Select job files to be migrated');
else
    vYjobs = varargin;
end
% Load jobs from files, pass others unchanged
LvYjobs = cfg_load_jobs(vYjobs);
%% Loop over jobs, migrate each
LvXjobs = cell(size(LvYjobs));
for cj = 1:numel(LvYjobs)
    if ~isempty(LvYjobs{cj})
        cLvXjob = local_migrate_fibertools_job(LvYjobs{cj});
        if ischar(vYjobs{cj}) % Job was loaded from file, save it to file
            [p n e] = fileparts(vYjobs{cj});
            LvXjobs{cj} = fullfile(p, [n '_vX' e]);
            jid = cfg_util('initjob', cLvXjob);
            cfg_util('savejob', jid, LvXjobs{cj});
            cfg_util('deljob', jid);
        else % Job was passed as variable
            LvXjobs{cj} = cLvXjob;
        end
    end
end
%% Return jobs
nout = min(nargout, numel(LvXjobs));
varargout = cell(1, nout);
[varargout{:}] = deal(LvXjobs{1:nout});

function vYjob = local_migrate_fibertools_job(vYjob)
% Loop over modules, try to replace old structures with new ones
for cm = 1:numel(vYjob)
    %% tracking.probabilistic - add choice to seed ROI and tracking area
    % old: seedmask and defmask contain mask number(s) or dependency
    % new: seedmask and defmask are choices - old inputs need to go into
    % field 'seedmask'/'defmask'
    try
        sm = vYjob{cm}.dtijobs.tracking.probabilistic.seedchoice.seedroi.seedmask;
        if ~isstruct(sm)
            vYjob{cm}.dtijobs.tracking.probabilistic.seedchoice.seedroi.seedmask = struct('seednumber', sm);
        end
    end
    try
        dm = vYjob{cm}.dtijobs.tracking.probabilistic.defroi.defarea.defmask;
        if ~isstruct(dm)
            vYjob{cm}.dtijobs.tracking.probabilistic.defroi.defarea.defmask = struct('defnumber', dm);
        end
    end
    %% tracking.mori - add choice to start and stop ROIs
    % old: startdef and stopdef contain mask number(s) or dependency
    % new: startdef and stopdef are choices - old inputs need to go into
    % field 'startmask'/'stopmask'
    try
        sm = vYjob{cm}.dtijobs.tracking.mori.start.startdef.startmask;
        if ~isstruct(sm)
            vYjob{cm}.dtijobs.tracking.mori.start.startdef.startmask = struct('startnumber', sm);
        end
    end
    try
        sm = vYjob{cm}.dtijobs.tracking.mori.stop.stopdef.stopmask;
        if ~isstruct(sm)
            vYjob{cm}.dtijobs.tracking.mori.stop.stopdef.stopmask = struct('stopnumber', sm);
        end
    end
    %% ftrop.selectbyROI
    % old: mask contains mask number(s) or dependency
    % new: mask is a choice - old inputs need to go into field 'masknumber'
    try
        sm = vYjob{cm}.dtijobs.ftrop.selectbyROI.roidef.mask;
        if ~isstruct(sm)
            vYjob{cm}.dtijobs.ftrop.selectbyROI.roidef.mask = struct('masknumber',sm);
        end
    end
    %% ftrop.eliminatebyROI
    % old: mask contains mask number(s) or dependency
    % new: mask is a choice - old inputs need to go into field 'masknumber'
    try
        sm = vYjob{cm}.dtijobs.ftrop.eliminatebyROI.roidef.mask;
        if ~isstruct(sm)
            vYjob{cm}.dtijobs.ftrop.eliminatebyROI.roidef.mask = struct('masknumber',sm);
        end
    end
    %% ftrop.visitMap
    % old: fibersubset is a subset number or dependency
    % new: fibersubset is a choice - old inputs need to go into field
    % 'fibernumber'
    try
        sm = vYjob{cm}.dtijobs.ftrop.visitMap.fibersubset;
        if ~isstruct(sm)
            vYjob{cm}.dtijobs.ftrop.visitMap.fibersubset = struct('fibernumber', sm);
        end
    end
    %% ftrop.visitMask
    % old: fibersubset is a subset number or dependency
    % new: fibersubset is a choice - old inputs need to go into field
    % 'fibernumber'
    try
        sm = vYjob{cm}.dtijobs.ftrop.visitMask.fibersubset;
        if ~isstruct(sm)
            vYjob{cm}.dtijobs.ftrop.visitMask.fibersubset = struct('fibernumber', sm);
        end
    end
    %% ftrop.endpointMask
    % old: fibersubset is a subset number or dependency
    % new: fibersubset is a choice - old inputs need to go into field
    % 'fibernumber'
    try
        sm = vYjob{cm}.dtijobs.ftrop.endpointMask.fibersubset;
        if ~isstruct(sm)
            vYjob{cm}.dtijobs.ftrop.endpointMask.fibersubset = struct('fibernumber', sm);
        end
    end
end

