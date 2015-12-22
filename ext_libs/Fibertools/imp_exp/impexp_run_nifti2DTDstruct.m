function varargout = impexp_run_nifti2dtdStruct(cmd, varargin)
% Template function to implement callbacks for an cfg_exbranch. The calling
% syntax is
% varargout = impexp_run_nifti2dtdStruct(cmd, varargin)
% where cmd is one of
% 'run'      - out = impexp_run_nifti2dtdStruct('run', job)
%              Run a job, and return its output argument
% 'vout'     - dep = impexp_run_nifti2dtdStruct('vout', job)
%              Examine a job structure with all leafs present and return an
%              array of cfg_dep objects.
% 'check'    - str = impexp_run_nifti2dtdStruct('check', subcmd, subjob)
%              Examine a part of a fully filled job structure. Return an empty
%              string if everything is ok, or a string describing the check
%              error. subcmd should be a string that identifies the part of
%              the configuration to be checked.
% 'defaults' - defval = impexp_run_nifti2dtdStruct('defaults', key)
%              Retrieve defaults value. key must be a sequence of dot
%              delimited field names into the internal def struct which is
%              kept in function local_def. An error is returned if no
%              matching field is found.
%              impexp_run_nifti2dtdStruct('defaults', key, newval)
%              Set the specified field in the internal def struct to a new
%              value.
% Application specific code needs to be inserted at the following places:
% 'run'      - main switch statement: code to compute the results, based on
%              a filled job
% 'vout'     - main switch statement: code to compute cfg_dep array, based
%              on a job structure that has all leafs, but not necessarily
%              any values filled in
% 'check'    - create and populate switch subcmd switchyard
% 'defaults' - modify initialisation of defaults in subfunction local_defs
% Callbacks can be constructed using anonymous function handles like this:
% 'run'      - @(job)impexp_run_nifti2dtdStruct('run', job)
% 'vout'     - @(job)impexp_run_nifti2dtdStruct('vout', job)
% 'check'    - @(job)impexp_run_nifti2dtdStruct('check', 'subcmd', job)
% 'defaults' - @(val)impexp_run_nifti2dtdStruct('defaults', 'defstr', val{:})
%              Note the list expansion val{:} - this is used to emulate a
%              varargin call in this function handle.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: impexp_run_nifti2DTDstruct.m,v 1.1 2012/12/04 13:48:16 reisertm Exp $

rev = '$Rev: 315 $'; %#ok

if ischar(cmd)
    switch lower(cmd)
        case 'run'
            job = local_getjob(varargin{1});
            % do computation, return results in variable out
            % create dtdStruct
            [dtd errstr] = nifti_to_DTD(job.srcimgsTensor, job.gdirTransformation);
            if ~isempty(errstr)
                error('impexp:dtdStruct',errstr);
            end
            switch char(fieldnames(job.outchoice))
                case 'outvar'
                    out.mskvar = dtd;
                case 'outmat'
                    [p n e] = fileparts(job.outchoice.outmat.fname);
                    out.mskmat = {fullfile(job.outchoice.outmat.outdir{1}, ...
                                          sprintf('%s.mat',n))};
                    [res errstr] = dtdstruct_write(dtd, out.mskmat{1});
                    if ~isempty(errstr)
                        error('impexp:dtdStruct',errstr);
                    end
            end
            if nargout > 0
                varargout{1} = out;
            end
        case 'vout'
            job = local_getjob(varargin{1});
            % initialise empty cfg_dep array
            dep = cfg_dep;
            % determine outputs, return cfg_dep array in variable dep
            switch char(fieldnames(job.outchoice))
                case 'outvar'
                    dep(1).sname   = 'dtdStruct variable';
                    dep.src_output = substruct('.','mskvar');
                    dep.tgt_spec   = cfg_findspec({{'strtype','e'}});
                case 'outmat'
                    dep(1).sname   = 'dtdStruct .mat file';
                    dep.src_output = substruct('.','mskmat');
                    dep.tgt_spec   = cfg_findspec({{'filter','mat', ...
                                        'strtype','e'}});
            end
            varargout{1} = dep;
        case 'check'
            if ischar(varargin{1})
                subcmd = lower(varargin{1});
                subjob = varargin{2};
                str = '';
                switch subcmd
                    % implement checks, return status string in variable str
                    otherwise
                        cfg_message('unknown:check', ...
                            'Unknown check subcmd ''%s''.', subcmd);
                end
                varargout{1} = str;
            else
                cfg_message('ischar:check', 'Subcmd must be a string.');
            end
        case 'defaults'
            if nargin == 2
                varargout{1} = local_defs(varargin{1});
            else
                local_defs(varargin{1:2});
            end
        otherwise
            cfg_message('unknown:cmd', 'Unknown command ''%s''.', cmd);
    end
else
    cfg_message('ischar:cmd', 'Cmd must be a string.');
end

function varargout = local_defs(defstr, defval)
persistent defs;
if isempty(defs)
    % initialise defaults
end
if ischar(defstr)
    % construct subscript reference struct from dot delimited tag string
    tags = textscan(defstr,'%s', 'delimiter','.');
    subs = struct('type','.','subs',tags{1}');
    try
        cdefval = subsref(local_def, subs);
    catch
        cdefval = [];
        cfg_message('defaults:noval', ...
            'No matching defaults value ''%s'' found.', defstr);
    end
    if nargin == 1
        varargout{1} = cdefval;
    else
        defs = subsasgn(defs, subs, defval);
    end
else
    cfg_message('ischar:defstr', 'Defaults key must be a string.');
end

function job = local_getjob(job)
if ~isstruct(job)
    cfg_message('isstruct:job', 'Job must be a struct.');
end