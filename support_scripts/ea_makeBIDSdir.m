function ea_makeBIDSdir(BIDSdir)
pipelines = {'coregistration','normalization','brainshift'};

%coregistration
for pipeline = 1:length(pipelines)
    bidsdir.mainDir = fullfile(BIDSdir,pipelines{pipeline});
    bidsdir.anatDir = fullfile(bidsdir.mainDir,'anat');
    bidsdir.checkregDir = fullfile(bidsdir.mainDir,'checkreg');
    bidsdir.transformationDir = fullfile(bidsdir.mainDir,'transformations');
    bidsdir.logDir = fullfile(bidsdir.mainDir,'log');
    dirnames = fieldnames(bidsdir);
    for i=1:length(dirnames)
        if ~exist(bidsdir.(dirnames{i}),'dir')
            mkdir(bidsdir.(dirnames{i}))
        end
    end
end
   
%reconstruction
bidsdir.other.recoDir = fullfile(BIDSdir,'reconstruction');
%preprocessing dir
bidsdir.other.preproDir = fullfile(BIDSdir,'preprocessing');
%log dir
bidsdir.other.logdir = fullfile(BIDSdir,'log');
%prefs dir
bidsdir.other.prefdir = fullfile(BIDSdir,'prefs');
new_dirname = fieldnames(bidsdir.other);
for i=1:length(new_dirname)
    if ~exist(bidsdir.other.(new_dirname{i}),'dir')
        mkdir(bidsdir.other.(new_dirname{i}))
        if strcmp(new_dirname{i},'preproDir')
            mkdir(bidsdir.other.(new_dirname{i}),'anat')
        end
    end
end


