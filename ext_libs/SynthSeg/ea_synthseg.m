function ea_synthseg(input, output)
% Wrapper to run SynthSeg
arguments
    input       % path to the nifti to segment
    output      % path to the segmented image
end

% Check Conda environment
condaenv = ea_conda_env('SynthSeg');
if ~condaenv.is_created
    ea_cprintf('*Comments', 'Initializing SynthSeg conda environment...\n')
    condaenv.create;
    ea_cprintf('*Comments', 'SynthSeg conda environment initialized.\n')
elseif ~condaenv.is_up_to_date
    ea_cprintf('*Comments', 'Updating SynthSeg conda environment...\n')
    condaenv.update;
    ea_cprintf('*Comments', 'SynthSeg conda environment initialized.\n')
end

% Run SynthSeg
synthseg_exe = fullfile(ea_getearoot, 'ext_libs', 'SynthSeg', 'mri_synthseg');

synthseg_cmd = {'python', ea_path_helper(synthseg_exe), ...
    '--i', ea_path_helper(input), ...
    '--o', ea_path_helper(output), ...
    '--parc --robust --threads -1 --noaddctab'};

status = condaenv.system(strjoin(synthseg_cmd, ' '));
if status ~= 0
    ea_error('SynthSeg failed!', showdlg=false, simpleStack=true);
end
