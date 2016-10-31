% mrInit.m
%
%      usage: scanParams = fixFramePeriods(scanParams,weirdFramePeriods,minFramePeriod,maxFramePeriod,tseriesDir)

%         by: julien besle, taken out of mrInit
%       date: 22/07/13
%    purpose: prompts user to fix apparently abnormal frame periods, called by mrInit and saveNewTSeries


function scanParams = fixFramePeriods(scanParams,weirdFramePeriods,minFramePeriod,maxFramePeriod,tseriesDir)

question{1} = 'Abnormal frame periods have been detected';
question{2} = 'Do you want to change the following values ?';
newFramePeriod = weirdFramePeriods;

for iScan = 1:length(weirdFramePeriods)
  if weirdFramePeriods(iScan)
    if weirdFramePeriods(iScan)<minFramePeriod
      newFramePeriod(iScan) = weirdFramePeriods(iScan)*1000;
    elseif weirdFramePeriods(iScan)>maxFramePeriod
      newFramePeriod(iScan) = weirdFramePeriods(iScan)/1000;
    end
    question{end+1} = sprintf('Change %f s into %f s in file %s',weirdFramePeriods(iScan),newFramePeriod(iScan),scanParams(iScan).fileName);
  end
end
yes = askuser(question);
if yes
  for iScan = 1:length(weirdFramePeriods)
    if newFramePeriod(iScan)
      filename = fullfile(tseriesDir, scanParams(iScan).fileName);
      hdr = mlrImageReadNiftiHeader(filename);
      %set time units to seconds (4th bit = 8)
      niftiSpaceUnit = bitand(hdr.xyzt_units,hex2dec('07')); 
      niftiTimeUnit = bitand(hdr.xyzt_units,hex2dec('38'));
      switch(niftiTimeUnit)
        case 8
          timeUnit = 'sec';
        case 16
          timeUnit = 'msec';
        case 32
          timeUnit = 'microsec';
      end
      % change units to seconds
      hdr.xyzt_units = bitor(bitand(niftiSpaceUnit,hex2dec('07')),bitand(8,hex2dec('38')));
      fprintf('Changing frame period from %f %s to %f sec in file %s\n',weirdFramePeriods(iScan),timeUnit,newFramePeriod(iScan),scanParams(iScan).fileName);
      % set the frameperiod
      scanParams(iScan).framePeriod = newFramePeriod(iScan);
      % change the nifti header
      hdr.pixdim(5) = newFramePeriod(iScan); 
      hdr.xyzt_units = 8+niftiSpaceUnit;
      % set in scanParams
      scanParams(iScan).niftiHdr = hdr;
      % and save back the header
      cbiWriteNiftiHeader(hdr,filename);
    end
  end
end
