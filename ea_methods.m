function ea_methods(options,parsestr,refs)
% function that dumps a methods text to patient directory
% options can either be a leadsuite options struct or a string with the
% patient directory.

if ~isstruct(options)
    directory=options; % direct supply of directory string
else
directory=[options.root,options.patientname,filesep];
end
if ~strcmp(directory(end),filesep) % make sure directory has a / in the end
    directory=[directory,filesep];
end

h=dbstack;
callingfunction=h(4).name;
expstr='\n\n';
expstr=[expstr,[datestr(datetime('now')),': ',callingfunction,'\n','--------------------------\n',...
    parsestr]];

if exist('refs','var') % add refs
    expstr=[expstr,'\n\nReferences:\n','--------------------------\n'];
    
    for r=1:length(refs)
        expstr=[expstr,[num2str(r),') ',refs{r},'\n']];
    end
end
expstr=[expstr,'\n\n***'];

fprintf(expstr);

metfile=fopen([directory,'ea_methods.txt'],'a');

fprintf(metfile,expstr);
fclose(metfile);


