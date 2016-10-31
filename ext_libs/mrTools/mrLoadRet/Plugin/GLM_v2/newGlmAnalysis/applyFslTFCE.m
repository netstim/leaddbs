% performs Threshold-Free Cluster Enhamcement on a 3-D array using flsmaths

function tfceData = applyFslTFCE(data, tempFilename,verbose)

if ieNotDefined('verbose')
  verbose = 1;
end
if ieNotDefined('tempFilename')
  %find temporary file extension based on FSL preference
  switch(getenv('FSLOUTPUTTYPE'))
    case 'NIFTI'
      tempFilename='temp.nii';
    case 'NIFTI_PAIR'
      tempFilename='temp.img';
    case ''
      mrWarnDlg('(fslApplyWarpSurfOFF) Environment variable FSLOUTPUTTYPE is not set');
      return;
    otherwise
      mrWarnDlg(sprintf('(fslApplyWarpOverlays) Environment variable FSLOUTPUTTYPE is set to an unsupported value (%s)',getenv('FSLOUTPUTTYPE')));
      return;
  end
end


dataSize = size(data);
if length(dataSize)<3 
  dataSize(3) = 1;
elseif length(dataSize)>4
  mrErrorDlg('(applyFslTFCE) Volume must be 3D or 4D');
end

newData = zeros([dataSize(1:3)+2 size(data,4)]);
newData(2:end-1,2:end-1,2:end-1,:) = data;

fslPath = mrGetPref('fslPath');
if strcmp(mrGetPref('fslPath'),'FSL not installed')
  mrErrorDlg('(applyFslTFCE) No path was provided for FSL. Please set MR preference ''fslPath'' by running mrSetPref(''fslPath'',''yourpath'')')
end

cbiWriteNifti(tempFilename, newData,[],[],[],[],verbose);
tfce_H = 2.0;
tfce_E = 0.5;
tfce_connectivity = 6;
try
  [s,w] = unix(sprintf('%s/fslmaths %s -tfce %.2f %.2f %d %s',fslPath,  tempFilename, tfce_H, tfce_E, tfce_connectivity, tempFilename));
  if s ~=- 0 % unix error
    disp('UNIX error message:')
    disp(w)
    disp('-------------------')
    return
  end
catch 
  disp('(applyFslTFCE) There was a problem running the TFCE unix command')
  disp(sprintf('unix error code: %d; %s', s, w))
  return
end
%read the TFCE values
tfceData=mlrImageReadNifti(tempFilename,[],[],[],verbose);

if all(tfceData(:)==0)
  oneTimeWarning('tfceOutputsZeros','(applyFslTFCE) There is a problem with fslmaths -tfce (it outputs only zeros). try using another version of FSL');
end

tfceData = tfceData(2:end-1,2:end-1,2:end-1,:);