function ea_run_HMS_Orchestra_medium(options)

% This is a function that runs code for each subject on a cluster. It needs
% to be adapted to suit your needs. Then, you can "Export code" using
% Lead-DBS and simply change the command ea_run to ea_run_cluster.

jobFile = [options.root, options.patientname, filesep, 'job_', ea_generate_uuid];
options.spmdir = spm('dir');
save(jobFile, 'options');

setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS', '1');

cmdstring=['cd ', options.earoot, ' &&' ...
           ' bsub -q medium -R "rusage[mem=70000]" -R "select[scratch2]" -W 119:0' ...
           ' -o ', jobFile,'.out' ...
           ' -e ', jobFile,'.err' ...
           ' matlab -singleCompThread -nodisplay -r "ea_run runcluster ', jobFile, '"'];

system(cmdstring);
