function [savedir]=ea_getfsstatesavedir(directory,vsname,fibersfile,seedfile,targetsfile,thresh,mode)
% helper function to save fiber states
% called from ea_cvshowvatdmri.m

[~,fibersname,~] = fileparts(fileparts(fibersfile));

if size(seedfile,2)==2
    sides='both';
else
    [~,sides,~]=fileparts(seedfile);
    sides=sides(length(mode)+2:end);
end

[~,parcname,~] = fileparts(targetsfile);

savedir = [fullfile(directory,'stimulations',vsname,'connvisfibers',fibersname,parcname,sides,thresh,mode),'/'];
