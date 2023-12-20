function tests = ea_explorer_statlist
    warning('off','all')
    tests = table('Size',[0,5],...
        'VariableTypes',{'categorical','string','categorical','cell','cell'},...
        'VariableNames',{'name','file','type','outcometype','compatibility'});
    row=0;

[pth]=fileparts(mfilename('fullpath')); % get dir
funs=dir(fullfile(pth,'ea_explorer_stats_*.m'));

for fun=1:length(funs)
entry=feval(ea_stripext(funs(fun).name));
tests.name(fun)=entry.name;
tests.file(fun)=entry.file;
tests.type(fun)=entry.type;
tests.outcometype{fun}=entry.outcometype;
tests.compatibility{fun}=entry.compatibility;
end

end