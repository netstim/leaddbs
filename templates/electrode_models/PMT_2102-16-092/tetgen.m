clear

smesh = ea_regexpdir(pwd, '.*\.smesh$');
for s=1:length(smesh)
    system(['"', mcpath('tetgen'), getexeext, '" -pq -aA -g "', smesh{s}, '"']);
end

ea_delete(ea_regexpdir(rootFolder, '.*\.edge$'));
ea_delete(ea_regexpdir(rootFolder, '.*\.mesh$'));
ea_delete(ea_regexpdir(rootFolder, '.*\.mtr$'));
