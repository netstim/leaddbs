function ea_gen_models(uipatdirs,typemap,regressor,fdmri,fis)
load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']); % here only needed to get a mask

odir=[fileparts(uipatdirs{1}),filesep,'predictive_model',filesep];
if ~exist(odir,'dir')
    mkdir(odir);
end

for pat=0
    opats=1:length(uipatdirs);
    if pat % zero model is leave-nothing-out model
        opats(pat)=[];
    end
    if ismember(fdmri.do,{'fmri','both'})
        feval(['ea_',typemap],fis.fmri(opats),regressor(opats),[odir,'fMRI_',typemap,'_leno_',num2str(pat),'.nii'],modeldata.mask);
    end
    if ismember(fdmri.do,{'dmri','both'})
        feval(['ea_',typemap],fis.dmri(opats),regressor(opats),[odir,'dMRI_',typemap,'_leno_',num2str(pat),'.nii'],modeldata.mask,'sk');
    end
end