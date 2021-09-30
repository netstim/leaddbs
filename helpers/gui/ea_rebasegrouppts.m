function ea_rebasegrouppts(~,~,handles)

newtarget=uigetdir('','Please select new base folder...');

Mfilebase=handles.groupdir_choosebox.String;
close(handles.leadfigure);

load([Mfilebase,'LEAD_groupanalysis.mat']);

for pt=1:length(M.patient.list)
    [pth,fn,ext]=fileparts(M.patient.list{pt});
    M.patient.list{pt}=fullfile(newtarget,fn);
end

save([Mfilebase,'LEAD_groupanalysis.mat'],'M','-v7.3');

lead group

