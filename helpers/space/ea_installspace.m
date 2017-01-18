function ea_installspace


disp(['Installing / Downloading space ',ea_getspace,'...']);
disp('This could take a while...');
downloadurl = 'http://www.lead-dbs.org/release/download.php';
    webopts=weboptions('Timeout',5);
    destination=[ea_space,'../data_download.zip'];
    try
        websave(destination,downloadurl,'id',ea_getspace,webopts);
    catch
        urlwrite([downloadurl,'?id=',ea_getspace],destination,'Timeout',5);
    end
disp('Download done. Will now continue building/unpacking space.');

space=ea_space;
rmdir(space,'s');
unzip(destination,space);
[~,spacename]=fileparts(fileparts(space));
movefile([space,spacename],[space(1:end-1),'_tmp']);
rmdir(space,'s');
movefile([space(1:end-1),'_tmp'],[space(1:end-1)]);
delete(destination);
ea_unpackspace;

