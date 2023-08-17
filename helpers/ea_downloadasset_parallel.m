function ea_downloadasset_parallel(downloadurl, assetname, destination, id, fsize)

%% doc
% EA_DOWNLOADASSET_PARALLEL utilizes the parallel toolbox to download a
% file and display its progress with a waitbar in parallel
%
%   input arguments:
%       downloadurl: str, url to be passed to websave
%       assetname: str, name of the asset to be downloaded, will be displayed in
%       the title of the waitbar figure
%       destination: str, location at which the file is to be saved
%       id: str, id of the file to be downloaded, passed to websave
%       fsize: double, file size in byte
%
%   Author: Johannes Achtzehn, 19.11.2020, johannes.achtzehn@charite.de

%% code
downloadaborted = false;        % this stays false if download is succesful, if waitbar is cancelled this will switch to true

% init 2 parallel workers, one for downloading, one for checking the
% progress
delete(gcp('nocreate'))         % delete any existing parallel workers
p = parpool(2);                 % create 2 workers
q = parallel.pool.DataQueue;    % create data queue so we can send data from the parallel workers

% create waitbar figure with cancel button that starts canceldownload
% function
wbhandle = waitbar(0, 'Starting download...', ...
        'Name', sprintf('Downloading %s', assetname), ...
        'CreateCancelBtn', @(obj, event) ea_canceldownload());

warn('off', 'MATLAB:parallel:loadsave:CannotSaveCorrectly');
f(1) = parfeval(p, @ea_startwebsave, 0, destination, downloadurl, id);      % start the download
f(2) = parfeval(p, @ea_checkprogress, 0, q, destination, fsize);            % start monitoring of file size

% after each catches updates from ea_checkprogress and updates the waitbar
afterEach(q, @(i) waitbar(i/100, wbhandle, ...
    sprintf('%.2f%% (%.2f/%.2f GB)', i, i/100*fsize*1e-9, fsize*1e-9)));

wait(f(1));     % wait for download to be completed

% if everything went well, downloadaborted is still false
if ~downloadaborted    
    waitbar(1, wbhandle, 'Finished downloading');
    ea_cleanup();
% if download was aborted via the cancel button, throw an error that is
% catched outside of this function
else
    error('Parallel download unsuccessful')
end

warn('on', 'MATLAB:parallel:loadsave:CannotSaveCorrectly');

% helper functions
    function ea_checkprogress(q, filename, targetFilesize) 

        currentFilesize = 0;    % init variable that tracks file size in bytes
        
        % wait for file to be created until checking its file size
        while isfile(filename) == 0
            pause(1)
        end
        
        % while current filesize si smaller than target file size, check
        % periodically size and send fractional size
        while currentFilesize < targetFilesize
            s = dir(filename);
            currentFilesize = s.bytes;
            i = (currentFilesize/targetFilesize) * 100;
            send(q, i);
            pause(0.5)
        end
    end 
        
    function ea_startwebsave(destination, downloadurl, id)
        webopts=weboptions('Timeout',Inf);
        websave(destination,downloadurl,'id',id, webopts);
    end

    function ea_cleanup()
        delete(gcp('nocreate'));    % shutdown workers
        % delete waitbar
        try 
            delete(findall(0,'type','figure','tag','TMWWaitbar'));
        catch
            warning('Tried to delete waitbar figure but found none.')
        end
    end

    function ea_canceldownload()
        delete(gcp('nocreate'));  
        downloadaborted = true;   
    end
end