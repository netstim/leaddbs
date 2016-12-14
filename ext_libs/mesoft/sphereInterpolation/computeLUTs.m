function computeLUTs

dirs = load('dwidirections.mat');

fns = fieldnames(dirs);
for k = 1:length(fns),
    bDir = getfield(dirs,fns{k});
    lu = sphereInterpolLUT([bDir -bDir]');
    save(sprintf('sinterp%istruct',size(bDir,2)),'-struct','lu');
   
end;