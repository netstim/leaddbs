function ea_run_Orchestra(options)

% This is a function that runs code for each subject on a cluster. It needs
% to be adapted to suit your needs. Then, you can "Export code" using
% Lead-DBS and simply change the command ea_run to ea_run_cluster.


jobID=ea_generate_guid;
save([options.root,options.patientname,filesep,jobID],'options')
cmdstring=['bsub -q short -W 10:0 -o ',[options.root,options.patientname,filesep,jobID],'.out -e ',[options.root,options.patientname,filesep,jobID],'.err matlab -nodisplay -r "ea_run runcluster ',[options.root,options.patientname,filesep,jobID],'"'];
system(cmdstring);
