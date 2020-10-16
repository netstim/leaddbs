function ea_methods(options,parsestr,refs)
% function that dumps a methods text to patient directory
% options can either be a leadsuite options struct or a string with the
% patient directory.

h=dbstack;
try
    callingfunction=h(2).name;
catch
    callingfunction='base';
end

expstr='\n\n';
expstr=[expstr,[datestr(datetime('now')),': ',callingfunction,'\n','--------------------------\n',...
    parsestr]];

if exist('refs','var') % add refs
    expstr=[expstr,'\n\nReferences:\n','--------------------------\n'];

    for r=1:length(refs)
        expstr=[expstr,[num2str(r),') ',refs{r},'\n']];
    end
end

expstr=[expstr,'\n***\n\n'];

prefs = ea_prefs;
if prefs.machine.methods_show
    try
        ea_methodsdisp({expstr});
    end
else
    fprintf(expstr);
end

if ~isempty(options)
    if ischar(options)
        options=ea_getptopts(options);
    end

    % Use try...catch since options may not have root and patientname field
    try
        methodfile=fopen([options.root,options.patientname,filesep,'ea_methods.txt'],'a');
        fprintf(methodfile,expstr);
        fclose(methodfile);
    end
end

