clear

smesh = ea_regexpdir(pwd, '.*\.smesh$');
for s=1:length(smesh)
    system(['"', ea_getExec(mcpath('tetgen')), '" -pq -aA -g "', smesh{s}, '"']);
end

ea_delete(ea_regexpdir(pwd, '.*\.edge$'));
ea_delete(ea_regexpdir(pwd, '.*\.mesh$'));
ea_delete(ea_regexpdir(pwd, '.*\.mtr$'));
