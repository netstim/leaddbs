clear
lead path % configure Lead-DBS path
% R-map Leave-one-out cross validation example script
allpts=1:100; % all patients

usemask='cortexcb'; % which areas of the brain to compare when comparing left-out patient's connectivity map with R-map generated from all other patients. This assumes 2x2x2 MNI-Space resolution. If using 1 mm res, could use _hd extension (see ea_getmask).
corrtype='Spearman'; % use which type of correlation metric - for fMRI could also do Pearson

for pt=allpts % iterate through patients - this loop will just give us a cell with entries that map to each patients structural connectivity map (seeding from the VTAs).
    
   patConnectivityFiles{pt}=['/path/to/my/folder/patient_',sprintf('%02.0f',pt),'_struc.nii']; % point to the files seeding from VTAs (could also be functional connectivity and/or could be generated using Lead connectome mapper).
    
end


Regressor=randn(length(allpts),1); % replace this regressor with your outcome variable, needs to match order of patConnectivityFiles.
ea_mkdir('Rmaps');


for pt=[0,allpts] % iterate through patients, leaving out one each time. In first iteration ("0"), we will leave nothing out and generate an R-map over all patients (usually denoted "R0map" in our nomenclature).
    
    otherpts=allpts;
    if pt>0 % in all but the first iteration, delete one patient from "otherpts".
        otherpts(otherpts==pt)=[]; % delete patient
    end
    %% 1. generate R-Maps
    
    
    ea_Rmap(patConnectivityFiles(otherpts),Regressor(otherpts),['Rmaps',filesep,'R_',sprintf('%02.0f',pt),'.nii'],ea_getmask(usemask),'',corrtype);
    
    if pt==0 % visualize main R-map (across all patients) in surfice
        ea_surficeoverlay(['Rmaps',filesep,'R_',sprintf('%02.0f',pt),'.nii']);
    end
    
    
    %% 2. now compare that R-map with the left-out patient's connectivity
    if pt>0
        thisPatConn=ea_load_nii(patConnectivityFiles{pt});
        thisRmap=ea_load_nii(['Rmaps',filesep,'R_',sprintf('%02.0f',pt),'.nii']);
        
        RegressorHat(pt)=corr(thisPatConn.img(ea_getmask(usemask),thisRmap.img(ea_getmask(usemask))),'rows','pairwise','type',corrtype); % estimate of how similar this patient's connectivity is to the "optimal" connectivity profile denoted by the R-map (that is based on all patients except this particular one).
    end
end

%% 3. show correlation between similarities to "optimal" connectivity and empirical improvement

h=ea_corrplot(Regressor(allpts), RegressorHat(allpts), 'noperm', {'LOOCV Crossvalidation','Empirical Regressor','Similarity to R-Map'});
saveas(h,'result.png');





