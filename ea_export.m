function ea_export(options)
% This function exports jobs created by the GUI of Lead-DBS. It is
% distributed within Lead-DBS toolbox (www.lead-dbs.org)
% __________________________________________________________________________________
% Copyright (C) 2016 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn



[fn,pth]=uiputfile('*.m','Specify location for new job file...','lead_job.m');
if ~fn % user pressed cancel
    return
end
exp=ea_gencode(options,'options');

fID=fopen([pth,fn],'w');

% export comments
fprintf(fID,'%s\n',['% - Lead-DBS Job created on ',date,':']);
fprintf(fID,'%s\n',['% --------------------------------------']);
fprintf(fID,'\n');
fprintf(fID,'\n');

for e=1:length(exp)
    fprintf(fID,'%s\n',exp{e});
end
fprintf(fID,'\n');
fprintf(fID,'%s\n','% Execute job:');
fprintf(fID,'%s\n','% ---------------------------------------');
fprintf(fID,'\n');
fprintf(fID,'%s\n',['ea_run(''run'',options);']);

edit([pth,fn]);