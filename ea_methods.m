function ea_methods(options,parsestr,refs)
% Show method text and dump it to patient directory if possible
%
% options can either be a struct or a string with the patient directory.

h = dbstack;
try
    callingfunction = h(2).name;
catch
    callingfunction = 'base';
end

expstr = '\n\n';
expstr = [expstr,[datestr(datetime('now')),': ',callingfunction,'\n','--------------------------\n',...
    parsestr]];

if exist('refs','var') % add refs
    expstr = [expstr,'\n\nReferences:\n','--------------------------\n'];

    for r=1:length(refs)
        expstr = [expstr,[num2str(r),') ',refs{r},'\n']];
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

if ~isempty(options) && isfield(options, 'subj')
    if ischar(options)
        options = ea_getptopts(options);
    end

    ea_mkdir(fileparts(options.subj.methodLog));
    methodfile = fopen(options.subj.methodLog, 'a');
    fprintf(methodfile, expstr);
    fclose(methodfile);
end
