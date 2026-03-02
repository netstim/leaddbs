function ea_get_ANN_env

env = ea_conda_env('SynthSeg');
env.force_create
if ~env.is_created
    ea_cprintf('*Comments', 'Installing SynthSeg conda environment for TensorFlow...\n');
    env.force_create;
    ea_cprintf('*Comments', 'Done.\n');
end

env.system('pip3 install matplotlib');
