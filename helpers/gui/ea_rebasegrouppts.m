function ea_rebasegrouppts(~,~,handles)

newtarget=uigetdir('','Please select new base folder...');

Mfilebase=handles.groupdir_choosebox.String;
close(handles.leadfigure);

[~,mid]=fileparts(fileparts(Mfilebase));
gfi=dir([Mfilebase,'dataset-*_analysis-*',mid,'.mat']);

load([Mfilebase,gfi(1).name]);

for pt=1:length(M.patient.list)
    [pth,fn,ext]=fileparts(M.patient.list{pt});
    M.patient.list{pt}=fullfile(newtarget,fn);
end

save([Mfilebase,gfi(1).name],'M','-v7.3');

lead group