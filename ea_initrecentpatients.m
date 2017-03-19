function ea_initrecentpatients(handles,patsub)

if ~exist('patsub','var')
    patsub='patients';
end

earoot=ea_getearoot;
try
    load([earoot,'common',filesep,'ea_recentpatients.mat']);
catch
    fullrpts={['No recent ',patsub,' found']};
end
save([earoot,'common',filesep,'ea_recentpatients.mat'],'fullrpts');
ea_updaterecentpatients(handles,patsub);