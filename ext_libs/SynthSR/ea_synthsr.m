function ea_synthsr(input, output, opts)
% Wrapper to run SynthSR
arguments
    input {mustBeTextScalar}
    output {mustBeTextScalar}
    opts.ct {mustBeNumericOrLogical} = false % For CT scans in Hounsfield scale
    opts.sharpening {mustBeNumericOrLogical} = true % Unsharp masking by default
    opts.flipping {mustBeNumericOrLogical} = true % Flipping augmentation at test time by default
    opts.lowfield {mustBeNumericOrLogical} = false % Use model for low-field scans
end

% Check Conda environment
condaenv = ea_conda_env('SynthSR');
if ~condaenv.is_created
    ea_cprintf('*Comments', 'Initializing SynthSR conda environment...\n')
    condaenv.create;
    ea_cprintf('*Comments', 'SynthSR conda environment initialized.\n')
elseif ~condaenv.is_up_to_date
    ea_cprintf('*Comments', 'Updating SynthSR conda environment...\n')
    condaenv.update;
    ea_cprintf('*Comments', 'SynthSR conda environment initialized.\n')
end

% Run SynthSR
synthsr_exe = fullfile(ea_getearoot, 'ext_libs', 'SynthSR', 'mri_synthsr');

synthsr_cmd = {'python', ea_path_helper(synthsr_exe), ...
    '--i', ea_path_helper(input), ...
    '--o', ea_path_helper(output)};

if isfield(opts, 'ct') && opts.ct
    synthsr_cmd = [synthsr_cmd, '--ct'];
end

if isfield(opts, 'sharpening') && ~opts.sharpening
    synthsr_cmd = [synthsr_cmd, '--disable_sharpening'];
end

if isfield(opts, 'flipping') && ~opts.flipping
    synthsr_cmd = [synthsr_cmd, '--disable_flipping'];
end

if isfield(opts, 'lowfield') && opts.lowfield
    synthsr_cmd = [synthsr_cmd, '--lowfield'];
end

synthsr_cmd = [synthsr_cmd, '--threads -1'];

status = condaenv.system(strjoin(synthsr_cmd, ' '));
if status ~= 0
    ea_error('SynthSR failed!', showdlg=false, simpleStack=true);
end
