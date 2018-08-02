function ea_run_Martinos_Launchpad_compile(options)

% This is a function that runs code for each subject on a cluster. It needs
% to be adapted to suit your needs. Then, you can "Export code" using
% Lead-DBS and simply change the command ea_run to ea_run_cluster.

jobFile = [options.root, options.patientname, filesep, 'job_', ea_generate_uuid];
options.spmdir = spm('dir');
save(jobFile, 'options');

setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS', '1');

cmdstring=['cd ', options.earoot, ' &&' ...
           ' pbsubmit -q highio' ...
           ' -O ',jobFile,'.out', ...
           ' -E ', jobFile, '.err', ...
           ' -c "', [ea_getearoot, 'cluster', filesep, 'run_ea_run_binary.sh'], ...
           ' /usr/pubsw/common/matlab/8.6',...
           ' runcluster ''', ea_path_helper(jobFile), '''"'];

if ~exist([ea_getearoot, 'cluster', filesep, 'ea_run_binary'],'file')
    mcc('-m', 'ea_run', '-o', 'ea_run_binary');
    movefile('ea_run_binary', ['cluster', filesep, 'ea_run_binary']);
    movefile('run_ea_run_binary.sh', ['cluster', filesep, 'run_ea_run_binary.sh']);
end

keyboard
system(cmdstring);
