function ea_initrecentpatients(handles,patsub)
if ~exist('patsub','var')
    patsub='patients';
end
try
    load([ea_prefsdir, filesep, 'ea_recent', patsub, '.mat']);
catch
    fullrpts={['No recent ',patsub,' found']};
end
save([ea_prefsdir, filesep, 'ea_recent', patsub, '.mat'],'fullrpts');
ea_updaterecentpatients(handles,patsub);