function ea_downloadasset_parallel(downloadurl, assetname, destination, id, fsize)

downloadaborted = false;
wbhandle = waitbar(0, 'Starting parallel worker...', ...
        'Name', sprintf('Downloading %s', assetname), ..., 
        'CreateCancelBtn', @(obj, event) ea_canceldownload());
    
delete(gcp('nocreate'))         % delete any existing parallel workers
p = parpool(2);                 % create 2 workers
q = parallel.pool.DataQueue;    % create data queue so we can send data from the parallel workers

f(1) = parfeval(p, @ea_startwebsave, 0, destination, downloadurl, id);                 % start the download
f(2) = parfeval(p, @ea_checkprogress, 0, q, destination, fsize);    % start monitoring of file size
afterEach(q, @(i) waitbar(i/100, wbhandle, ...
    sprintf('%.2f%% (%.2f/%.2f GB)', i, i/100*fsize*1e-9, fsize*1e-9)));

wait(f(1));

if ~downloadaborted
    waitbar(1, wbhandle, 'Finished downloading');  
    ea_cleanup();
else
    error('Parallel download unsuccessful')
end

    function ea_checkprogress(q, filename, targetFilesize) 

        currentFilesize = 0;
        % wait for file to be created
        while isfile(filename) == 0
            pause(1)
        end

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