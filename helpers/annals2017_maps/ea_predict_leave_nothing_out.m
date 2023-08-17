function ea_predict_leave_nothing_out(uipatdirs,stimname,typemap,regressor,fdmri,dowhat)
if size(regressor,1)<size(regressor,2)
    regressor=regressor';
end

if ~exist('fdmri','var')
    fdmri.do='both'; % can be fmri, dmri or both
end

for pat=1:length(uipatdirs)
    [~, subPrefix] = fileparts([uipatdir{pt}, '_']);
    fConnName = ea_getConnLabel(fdmri.fmriconnectome);
    dConnName = ea_getConnLabel(fdmri.dmriconnectome);
    fis.fmri{pat}=fullfile(uipatdirs{pat},'stimulations',ea_nt(0),stimname,[subPrefix, 'sim-binary_model-simbio_seed-fMRI_conn-', fConnName, '_desc-AvgRFz_funcmap.nii']);
    fis.dmri{pat}=fullfile(uipatdirs{pat},'stimulations',ea_nt(0),stimname,[subPrefix, 'sim-binary_model-simbio_seed-dMRI_conn-', dConnName, '_strucmap.nii']);
end
if ~exist('dowhat','var')
    dowhat='both';
end
if ismember(dowhat,{'both','genmodels'});
    ea_gen_models(uipatdirs,typemap,regressor,fdmri,fis);
end
if ismember(dowhat,{'both','predictoutcome'});
    ea_predict_outcomes(uipatdirs,stimname,typemap,regressor,fdmri,fis);
end


function ea_predict_outcomes(uipatdirs,stimname,typemap,regressor,fdmri,fis)
load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']); % here only needed to get a mask
odir=[fileparts(uipatdirs{1}),filesep,'predictive_model',filesep];

for pat=1:length(uipatdirs)
    if ismember(fdmri.do,{'fmri','both'})
        model=ea_load_nii([odir,'fMRI_',typemap,'_leno_',num2str(0),'.nii']);
        patcon=ea_load_nii(fis.fmri{pat});
        pred.fmri(pat)=corr(model.img(modeldata.mask),patcon.img(modeldata.mask),'rows','pairwise','type','spearman');
    end
    if ismember(fdmri.do,{'dmri','both'})
        model=ea_load_nii([odir,'dMRI_',typemap,'_leno_',num2str(0),'.nii']);
        patcon=ea_load_nii(fis.dmri{pat});
        pred.dmri(pat)=corr(model.img(modeldata.mask),patcon.img(modeldata.mask),'rows','pairwise','type','spearman');
    end
end

if ismember(fdmri.do,{'fmri','both'})
    [R_upd,p_upd,R,p,h]=ea_corrplot_gen([regressor,pred.fmri'],'Corr',{'Empirical Improvement';'Predicted Improvement'},[],[],'permutation_spearman');
    axis square
    ea_screenshot([odir,'leno_prediction_fmri_',typemap,'.png']);
end
if ismember(fdmri.do,{'dmri','both'})
    [R_upd,p_upd,R,p,h]=ea_corrplot_gen([regressor,pred.dmri'],'Corr',{'Empirical Improvement';'Predicted Improvement'},[],[],'permutation_spearman');
    axis square
    ea_screenshot([odir,'leno_prediction_dmri_',typemap,'.png']);
end
save([odir,'leno_prediction'],'pred');

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

