function ea_run_Martinos_Launchpad_compile(options)

% This is a function that runs code for each subject on a cluster. It needs
% to be adapted to suit your needs. Then, you can "Export code" using
% Lead-DBS and simply change the command ea_run to ea_run_cluster.


jobID=ea_generate_guid;
options.spmdir=spm('dir');
save([options.root,options.patientname,filesep,'job_',jobID],'options')
setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS','1')

cmdstring=['pbsubmit -q highio -O ',[options.root,options.patientname,filesep,'job_',jobID],'.out',...
    ' -E ',[options.root,options.patientname,filesep,'job_',jobID],'.err',...
    ' -c "',[ea_getearoot,filesep,'cluster',filesep,'ea_run_binary.',computer('arch')],'',...
    ' ','/usr/pubsw/common/matlab/8.6',...
    ' runcluster ',[options.root,options.patientname,filesep,'job_',jobID],'"'];

if ~exist([ea_getearoot,'cluster',filesep,'ea_run_binary.',computer('arch')],'file');
mcc('-m','ea_run','-o',['ea_run_binary']);
movefile('ea_run_binary',['cluster',filesep,'ea_run_binary.',computer('arch')]);
delete('run_ea_run_binary.sh');
end
keyboard
system(cmdstring);
