function ea_get_ANN_env

env = ea_conda_env('SynthSeg');
if ~env.is_created
    ea_cprintf('*Comments', 'Installing SynthSeg conda environment for TensorFlow...\n');
    env.force_create;
    ea_cprintf('*Comments', 'Done.\n');
end

env.system('pip3 install matplotlib');
env.system('pip3 install seaborn');

% Set python path
binPath = getenv('PATH');
if isunix
    pythonPath = [env.path, filesep, 'bin'];
    setenv('PATH', [pythonPath, ':', binPath]);
else
    pythonPath = [env.path,';',env.path,filesep,'Scripts'];
    setenv('PATH', [pythonPath, ';', binPath]);
end