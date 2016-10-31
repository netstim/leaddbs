function [view] = saveTSeries(view,tseries,scanNum,scanParams,hdr,append)
%
% view = saveTSeries(view,tseries,scanNum,[scanParams],[hdr],[append])
%
% Note, this function is used to OVERWRITE a tSeries. If you want
% to save a new tSeries, you should use saveNewTSeries.
%
% tseries: tseries array (x,y,z,t)
% scanNum: the scan number to write to
% scanParams: structure with the following fields: fileName, description,
%    junkFrames, nFrames. The other scanParams fields are set automatically
%    from the nifti header. Default: Use the old scanParams.
%    If scanParams.fileName is not set, uses a default filename
%    'yymmdd-HHMMSS' and chooses extension (.img or .nii) based on
%    'niftiFileExtension' preference.  
% append: (1 or 0) Append the passed in time series to the existing time
%    series. Default: 0.
% hdr: template for nifti header. The header is always passed through
%    cbiCreateNiftiHeader to ensure consistency with the data. Default: [].

% check to make sure the scan exists
% get the old one
oldScanParams = viewGet(view,'scanParams',scanNum);
if isempty(oldScanParams)
    disp(sprintf('saveTSeries: Scan #%i does not exist',scanNum));
end

% get defaults
if ieNotDefined('scanParams'),scanParams = oldScanParams;end
if ieNotDefined('hdr'),hdr = [];end
if ieNotDefined('append'),append = 0;end

% default filename
ext = mrGetPref('niftiFileExtension');
if isempty(scanParams.fileName)
    scanParams.fileName = ['tseries-',datestr(now,'yymmdd-HHMMSS'),ext];
end
filename = scanParams.fileName;

% get the directory for the tseries
tseriesdir = viewGet(view,'tseriesdir');
path = fullfile(tseriesdir,scanParams.fileName);

if ~append
    % Save tseries
    [byteswritten,hdr] = cbiWriteNifti(path,tseries,hdr);
    scanParams.niftiHdr = hdr;
else
    % append time series to what is already there.
    % first get how many total frames we will have
    oldNFrames = viewGet(view,'nFrames',scanNum);
    newNFrames = oldNFrames+size(tseries,4);
    % get the old header if we don't have one passed in
    if isempty(hdr)
        hdr = oldScanParams.niftiHdr;
    end
    % and set how many total frames we will have.
    hdr.dim(5) = newNFrames;
    scanParams.totalFrames = newNFrames;
    scanParams.nFrames = newNFrames;
    % now write out that header
    hdr = cbiWriteNiftiHeader(hdr,path);
    scanParams.niftiHdr = hdr;
    % now write out the new data
    tempFilename = fullfile(tseriesdir,'___saveTSeriesAppendTemp___.img');
    cbiWriteNifti(tempFilename,tseries,hdr);
    % now append that file to the end of the old one
    % make sure to get the filename with .img appended
    [oldpath,oldname,oldext] = fileparts(path);
    oldImgFilename = fullfile(oldpath,sprintf('%s%s',oldname,oldext));
    % now append the files together
    fNew = fopen(tempFilename);
    if (fNew == -1)
        error(sprintf('saveTSeries: Could not open temporary file %s',tempFilename));
    end
    fOld = fopen(oldImgFilename,'a');
    if (fOld == -1)
        error(sprintf('saveTSeries: Could not open file %s',oldImgFilename));
    end
    % now cycle through and append blocks to the old file
    BLOCKSIZE = 512;
    [datablock count] = fread(fNew,BLOCKSIZE,hdr.matlab_datatype);
    while (count > 0)
        fwrite(fOld,datablock,hdr.matlab_datatype);
        [datablock count] = fread(fNew,BLOCKSIZE,hdr.matlab_datatype);
    end
    % and close
    fclose(fNew);fclose(fOld);
    % and delete temporary file
    delete(tempFilename);
    [tempPath tempFilename] = fileparts(tempFilename);
    delete(fullfile(tempPath,sprintf('%s.hdr',tempFilename)));
end

% Save scan params
view = viewSet(view,'updateScan',scanParams,scanNum);
