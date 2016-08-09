function ea_run_HMS_Orchestra(options)

% This is a function that runs code for each subject on a cluster. It needs
% to be adapted to suit your needs. Then, you can "Export code" using
% Lead-DBS and simply change the command ea_run to ea_run_cluster.


jobID=ea_generate_guid;
options.spmdir=spm('dir');
save([options.root,options.patientname,filesep,'job_',jobID],'options')
cmdstring=['bsub -q short -W 10:0 -o ',[options.root,options.patientname,filesep,'job_',jobID],'.out -e ',[options.root,options.patientname,filesep,'job_',jobID],'.err matlab -nodisplay -r "ea_run runcluster ',[options.root,options.patientname,filesep,'job_',jobID],'"'];
system(cmdstring);
