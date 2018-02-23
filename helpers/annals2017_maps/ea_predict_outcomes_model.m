function ea_predict_outcomes_model(uipatdirs,stimname,map,regressor,fdmri,fis)
load([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'modeldata.mat']); % here only needed to get a mask
odir=[fileparts(uipatdirs{1}),filesep,'predictive_model',filesep];

for pat=1:length(uipatdirs)
    if ismember(fdmri.do,{'fmri','both'})
        model=ea_load_nii(map.fmri);
        patcon=ea_load_nii(fis.fmri{pat});
        pred.fmri(pat)=corr(model.img(modeldata.mask),patcon.img(modeldata.mask),'rows','pairwise','type','spearman');
    end
    if ismember(fdmri.do,{'dmri','both'})
        model=ea_load_nii(map.dmri);
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

