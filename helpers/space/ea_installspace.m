function ea_installspace


disp(['Installing / Downloading space ',ea_getspace,'...']);
disp('This could take a while...');
downloadurl = 'http://www.lead-dbs.org/release/download.php';
    webopts=weboptions('Timeout',5);
    destination=[ea_space,'data.zip'];
    try
        websave(destination,downloadurl,'id',ea_getspace,webopts);
    catch
        urlwrite([downloadurl,'?id=',ea_getspace],destination,'Timeout',5);
    end
disp('Download done. Will now continue building/unpacking space.');
unzip(destination,ea_space);
delete(destination);

ea_unpackspace;

