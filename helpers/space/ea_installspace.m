function ea_installspace

disp(['Installing / Downloading space ',ea_getspace,'...']);
disp('This could take a while...');
downloadurl = 'https://www.lead-dbs.org/release/download.php';
    destination=[ea_space,'../data_download.zip'];
    try
        webopts=weboptions('Timeout',Inf);
        websave(destination,downloadurl,'id',ea_getspace,webopts);
    catch
        urlwrite([downloadurl,'?id=',ea_getspace],destination,'Timeout',Inf);
    end
disp('Download done. Will now continue building/unpacking space.');

unzip(destination,fileparts(fileparts(ea_space)));
% delete 'need_install' in user environment
% keep 'need_install' in dev environment
if ~exist([ea_getearoot,'.git'],'dir')
    ea_delete([ea_space,'need_install']);
end
delete(destination);
ea_unpackspace;
