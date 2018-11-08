function ea_writeMvat(M,nii,options,data,datalabels)
hshid=ea_datahash(M.ui.listselect);
percthresh=0.1;
if ~exist([options.root,options.patientname,filesep,'statvat_results',filesep,'models'],'dir')
    mkdir([options.root,options.patientname,filesep,'statvat_results',filesep,'models'])
end
nii.dt(1)=16;
for f=1:length(data)
    nii.img(:)=data{f}(:);
    nii.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'models',filesep,'statvat_',M.clinical.labels{M.ui.clinicallist},'_',datalabels{f},'_',hshid,'.nii'];
    ea_write_nii(nii);
    if ismember('N',datalabels) && ~strcmp(datalabels{f},'N')
        [~,ix]=ismember('N',datalabels);
        N=data{ix};
        nii.img(N<(max(N)*percthresh))=nan;
        nii.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'models',filesep,'statvat_',M.clinical.labels{M.ui.clinicallist},'_',datalabels{f},'_nthresh_',hshid,'.nii'];
        ea_write_nii(nii);
    end
end