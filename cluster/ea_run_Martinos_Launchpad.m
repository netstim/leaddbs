function ea_run_Orchestra(options)

% This is a function that runs code for each subject on a cluster. It needs
% to be adapted to suit your needs. Then, you can "Export code" using
% Lead-DBS and simply change the command ea_run to ea_run_cluster.


jobID=ea_generate_guid;
save(jobID,'options')
cmdstring=['bsub -q short -W 10:0 matlab -nodisplay -r "ea_run runcluster ',[pwd,filesep,jobID],'"'];
system(cmdstring);
