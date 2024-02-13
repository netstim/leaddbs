function ea_synthstrip(input, output, opts)
% Wrapper to run SynthStrip
arguments
    input {mustBeTextScalar}
    output {mustBeTextScalar}
    opts.mask {mustBeTextScalar}
    opts.csf {mustBeNumericOrLogical} = true
end

% Check Conda environment
condaenv = ea_conda_env('SynthStrip');
if ~condaenv.is_created
    ea_cprintf('CmdWinWarnings', 'Initializing SynthStrip conda environment...\n')
    condaenv.create;
    ea_cprintf('CmdWinWarnings', 'SynthStrip conda environment initialized.\n')
elseif ~condaenv.is_up_to_date
    ea_cprintf('CmdWinWarnings', 'Updating SynthStrip conda environment...\n')
    condaenv.update;
    ea_cprintf('CmdWinWarnings', 'SynthStrip conda environment initialized.\n')
end

% Run SynthStrip
synthstrip_exe = fullfile(ea_getearoot, 'ext_libs', 'SynthStrip', 'mri_synthstrip');

synthstrip_cmd = {'python', ea_path_helper(synthstrip_exe), ...
    '--i', ea_path_helper(input), ...
    '--o', ea_path_helper(output)};

if isfield(opts, 'mask')
    synthstrip_cmd = [synthstrip_cmd, '--mask', ea_path_helper(opts.mask)];
end

if isfield(opts, 'csf') && ~opts.csf
    synthstrip_cmd = [synthstrip_cmd, '--no-csf'];
end

status = condaenv.system(strjoin(synthstrip_cmd, ' '));
if status ~= 0
    ea_error('SynthStrip failed!', showdlg=false, simpleStack=true);
end
