function ea_writeMvat(M,nii,options,data,datalabels)
hshid=ea_datahash(M.ui.listselect);

for f=1:length(data)
nii.img(:)=data{f}(:);
nii.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_',M.clinical.labels{M.ui.clinicallist},'_',datalabels{f},'_',hshid,'.nii'];
ea_write_nii(nii);
end