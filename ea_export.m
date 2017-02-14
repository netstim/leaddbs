function ea_export(options,clusterfunctionname)
% This function exports jobs created by the GUI of Lead-DBS. It is
% distributed within Lead-DBS toolbox (www.lead-dbs.org)
% __________________________________________________________________________________
% Copyright (C) 2016 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn



[fn,pth]=uiputfile('*.m','Specify location for new job file...','lead_job.m');
if ~fn % user pressed cancel
    return
end
try
    options=rmfield(options,'root');
    options=rmfield(options,'patientname');
end
exp=ea_gencode(options,'options');

fID=fopen([pth,fn],'w');

% export comments

fprintf(fID,'%s\n',['function ',fn(1:end-2)]);

fprintf(fID,'%s\n',['% - Lead-DBS Job created on ',date,':']);
fprintf(fID,'%s\n',['% --------------------------------------']);
fprintf(fID,'\n');
fprintf(fID,'\n');

fprintf(fID,'\n');
fprintf(fID,'%s\n','% Execute job:');
fprintf(fID,'%s\n','% ---------------------------------------');
fprintf(fID,'\n');

fprintf(fID,'%s\n','options=getoptslocal;');
fprintf(fID,'\n');

fprintf(fID,'%s\n',['addpath(genpath(''',options.earoot(1:end-1),'''));']);
fprintf(fID,'%s\n',['addpath(''',fileparts(which('spm')),''');']);
fprintf(fID,'\n');
if exist('clusterfunctionname','var') % submit to cluster instead of directly running
    fprintf(fID,'%s\n',['% options.uipatdirs=ea_checknoerrorfolders(options.uipatdirs); % uncomment this line to only apply to folders that ran into errors the last time.']);
    fprintf(fID,'%s\n',['% options.uipatdirs=ea_checknofilefolders(options.uipatdirs,filename); % uncomment this line to only apply to folders that MISS a certain file.']);
    fprintf(fID,'%s\n',['% options.uipatdirs=ea_checkfilefolders(options.uipatdirs,filename); % uncomment this line to only apply to folders that CONTAIN a certain file.']);
end
fprintf(fID,'\n');

if ~isempty(options.uipatdirs)
    

    fprintf(fID,'%s\n','allpatdirs=options.uipatdirs;');    
    fprintf(fID,'%s\n','for pat=1:length(allpatdirs)');
    fprintf(fID,'%s\n','% set subject specific options:');
    fprintf(fID,'%s\n','options.root=[fileparts(allpatdirs{pat}),filesep];');
    fprintf(fID,'%s\n','[~,thispatdir]=fileparts(allpatdirs{pat});');
    fprintf(fID,'%s\n','options.patientname=thispatdir;');
    fprintf(fID,'%s\n','options.uipatdirs=allpatdirs{pat};');
    if exist('clusterfunctionname','var') % submit to cluster instead of directly running
        fprintf(fID,'%s\n',['clusterfunctionname=''',clusterfunctionname,''';']);
        fprintf(fID,'%s\n','feval(eval([''@'',clusterfunctionname]),options)');
    else
        fprintf(fID,'%s\n',['ea_run(''run'',options);']);
    end
    fprintf(fID,'%s\n','end');
    
else % no patient mode (also e.g. connectome mapper)
    if exist('clusterfunctionname','var') % submit to cluster instead of directly running
        fprintf(fID,'%s\n',['clusterfunctionname=''',clusterfunctionname,''';']);
        fprintf(fID,'%s\n','feval(eval([''@'',clusterfunctionname]),options)');
    else
        fprintf(fID,'%s\n',['ea_run(''run'',options);']);
    end
    
end
fprintf(fID,'\n');
fprintf(fID,'\n');
fprintf(fID,'\n');
fprintf(fID,'\n');



fprintf(fID,'%s\n',['function options=getoptslocal']);


for e=1:length(exp)
    fprintf(fID,'%s\n',exp{e});
end

edit([pth,fn]);
