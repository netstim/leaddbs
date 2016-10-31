function combineTransformOverlays(thisView)
% combineTransformOverlays(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%   combines (masked) Overlays according to matlab or custom operators  in current view and current analysis
%
% jb 07/07/2010
%
% $Id$

inputOutputTypeMenu = {'3D Array','4D Array','Scalar','Structure'};
combinationModeMenu = {'Apply function to all overlays','Apply function to each overlay','Recursively apply to overlay pairs'};

%default params
%get names of combine Functions in combineFunctions directory
functionsDirectory = [fileparts(which('combineTransformOverlays')) '/combineTransformOverlayFunctions/'];
combineFunctionFiles =  dir([functionsDirectory '*.m']);
for iFile=1:length(combineFunctionFiles)
  combineFunctions{iFile} = stripext(combineFunctionFiles(iFile).name);
end

combineFunctionsMenu = [{'User Defined'} combineFunctions];
params.customCombineFunction = '';%(''@(x)max(0,-norminv(x))';
params.nOutputOverlays = 1;
params.additionalArrayArgs = '';
params.additionalArgs = '';
params.clip = 0;
params.alphaClip = 0;
params.passView = 0;
params.baseSpace = 0;
baseSpaceInterpMenu = {'Same as display','nearest','linear','spline','cubic'};
params.exportToNewGroup = 0;
params.outputName = '';
if viewGet(thisView,'basetype')~=1
  baseSpaceOption='enable=0';
else
  baseSpaceOption='enable=1';
end

askForParams = 1;
while askForParams
  params = {...
   {'combineFunction',combineFunctionsMenu,'type=popupmenu','name of the function to apply. This is a list of existing combine functions in the combineFunctions directory. To use another function, select ''User Defined'' and type the function name below'},...
   {'customCombineFunction',params.customCombineFunction,'name of the function to apply. You can use any type of matlab function (including custom) that accepts either scalars or multidimensional arrays. Any string beginning with an @ will be considered an anonymous function and shoulde be of the form @(x)func(x), @(x,y)func(x,y) ..., where the number of variables equals the number of overlay inputs and additional arguments. '},...
   {'inputOutputType',inputOutputTypeMenu,'type=popupmenu','Type of arguments accepted by the combination function. ''3D Array'' will pass each input overlay as a 3D array. ''Scalar'' will apply the function to each element of the input overlay(s). ''4D Array'' wil concatenate overlays on the 4th dimension and pass the 4D array as a single argument to the function. 3D and 4D array are faster that but not all functions accept multidimensional arrays as inputs. Use ''4D Array'' for functions that operate on one dimension of an array (e.g. mean) and specify the dimension as an additional scalar argument (usually 4). Choose ''Structure'' to pass the whole overlay structure'},...
   {'combinationMode',combinationModeMenu,'type=popupmenu', 'How the selected overlays are input ot the combineFunction. If ''all'', all the selected overlays are given as input at once (the number of inputs expected by the function must match the number of selected overlays). If ''each'', the combine function is run separately for each overlay and must accept only one input overlay). If ''pair'', the combineFunction is run on pairs of consecutive selected overlays and must accept two input overlays.'},...
   {'nOutputOverlays',params.nOutputOverlays,'incdec=[-1 1]','round=1','minmax=[0 Inf]','Number of outputs of the combineFunction'},...
   {'additionalArrayArgs',params.additionalArrayArgs,'constant arguments for functions that accept them. Arguments must be separated by commas. for Array input/output type, each argument will be repeated in a matrix of same dimensions of the overlay '},...
   {'additionalArgs',params.additionalArgs,'constant scalar arguments for functions that take both arrays and scalars. These arguments will be input at the end of each function call. They must be separated by commas. If using vector argument, do not use commas inside square brackets. '},...
   {'passView',params.passView,'type=checkbox','Check this if the function requires the current mrLoadRet view'},...
   {'clip',params.clip,'type=checkbox','Mask overlays according to clip values'},...
   {'alphaClip',params.alphaClip,'type=checkbox','Mask overlays according to alpha overlay clip values'},...
   {'baseSpace',params.baseSpace,'type=checkbox',baseSpaceOption,'Transforms overlays into the current base volume before applying the transform/combine function, and back into overlay space afterwards. Only implemented for flat maps (all cortical depths are used).'},...
   {'baseSpaceInterp',baseSpaceInterpMenu,'type=popupmenu','contingent=baseSpace','Type of base space interpolation '},...
   {'exportToNewGroup',params.exportToNewGroup,'type=checkbox','contingent=baseSpace','Exports results in base sapce to new group, scan and analysis. Warning: for flat maps, the data is exported to a volume in an arbitrary space. ROIs and overlays defined outside this new group will not be in register.'},...
   {'outputName',params.outputName,'radical of the output overlay names'},...
   {'printHelp',0,'type=pushbutton','callback',@printHelp,'passParams=1','buttonString=Print combineFunction Help','Prints combination function help in command window'},...
          };
  params = mrParamsDialog(params, 'Overlay Combination parameters');
  % Abort if params empty
  if ieNotDefined('params'),return,end
  if strcmp(params.combineFunction,'User Defined')
    params.combineFunction = params.customCombineFunction;
  end 
  if isempty(params.baseSpace)
    params.baseSpace=0;
  end
  if isempty(params.exportToNewGroup)
    params.exportToNewGroup=0;
  end
  if strcmp(params.baseSpaceInterp,'Same as display')
    baseSpaceInterp='';
  else
    baseSpaceInterp=params.baseSpaceInterp;
  end
  
  inputOutputTypeMenu = putOnTopOfList(params.inputOutputType,inputOutputTypeMenu);
  combinationModeMenu = putOnTopOfList(params.combinationMode,combinationModeMenu);
  combineFunctionsMenu = putOnTopOfList(params.combineFunction,combineFunctionsMenu);
  baseSpaceInterpMenu = putOnTopOfList(params.baseSpaceInterp,baseSpaceInterpMenu);
  
  if strcmp(params.combinationMode,'Recursively apply to overlay pairs') && params.combineFunction(1)=='@'
    mrWarnDlg('(combineTransformOverlays) Anonymous functions cannot be applied recursively.');
  elseif isempty(params.combineFunction) || (strcmp(params.combineFunction,'User Defined') && isempty(params.customCombineFunction))
    mrWarnDlg('(combineTransformOverlays) Please choose a combination/transformation function.');
  elseif (params.clip || params.alphaClip) && params.baseSpace
    mrWarnDlg('(combineTransformOverlays) Base space conversion is not yet compatible with using (alpha) masking.');
  %elseif
    %other controls here
  else
    askForParams = 0;
    overlayList = selectInList(thisView,'overlays');
    if isempty(overlayList)
       askForParams = 1;
    end
  end
end
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;

%get the overlay data
nScans = viewGet(thisView,'nScans');   
overlayData = viewGet(thisView,'overlays');
overlayData = overlayData(overlayList);
if params.baseSpace
  base2scan = viewGet(thisView,'base2scan');
  if any(any(abs(base2scan - eye(4))>1e-6)) || viewGet(thisView,'basetype')>0 %check if we're in the scan space
    baseCoordsMap=cell(nScans,1);
    %if not, transform the overlay to the base space
    for iScan = 1:nScans
      %here could probably put all overlays of a single scan in a 4D array, but the would have to put it back 
      % into the overlays structure array if inputOutputType is 'structure' 
      for iOverlay = 1:length(overlayData) 
        [overlayData(iOverlay).data{iScan}, voxelSize, baseCoordsMap{iScan}] = getBaseSpaceOverlay(thisView, overlayData(iOverlay).data{iScan},[],[],baseSpaceInterp);
      end
    end
  end
end

% mask the overlay data if needed
%   Rk: to implement conversion to base space with masking/clip masking, it would be better
%   to get both the overlays and masks at the same time and use maskoverlay
%   with the boxInfo option
if params.clip
   mask = maskOverlay(thisView,overlayList);
   for iScan = 1:length(mask)
      for iOverlay = 1:length(overlayData)
        if ~isempty(overlayData(iOverlay).data{iScan})
          overlayData(iOverlay).data{iScan}(~mask{iScan}(:,:,:,iOverlay))=NaN;
        end
      end
   end
end
if params.alphaClip
  alphaOverlayNum = zeros(1,length(overlayData));
  for iOverlay = 1:length(overlayData)
    alphaOverlayNum(iOverlay) = viewGet(thisView,'overlaynum',overlayData(iOverlay).alphaOverlay);
  end
  mask = maskOverlay(thisView,alphaOverlayNum);
  for iScan = 1:length(mask)
    for iOverlay = 1:length(overlayData)
      if alphaOverlayNum(iOverlay) &&  ~isempty(overlayData(iOverlay).data{iScan})
        overlayData(iOverlay).data{iScan}(~mask{iScan}(:,:,:,iOverlay))=NaN;
      end
    end
  end
end

%overlay names
for iOverlay = 1:length(overlayList)
  overlayNames{iOverlay} = overlayData(iOverlay).name;
end

%reformat input data
switch(params.inputOutputType)
  case 'Structure'
    overlayData = num2cell(overlayData);
  case {'3D Array','4D Array','Scalar'}
    newOverlayData = cell(nScans,length(overlayList));
    for iOverlay = 1:length(overlayList)
      for iScan = 1:nScans
        newOverlayData{iScan,iOverlay} = overlayData(iOverlay).data{iScan};
      end
    end
    overlayData = newOverlayData;
end

%parse additional array inputs
additionalArrayArgs = parseArguments(params.additionalArrayArgs,',');
if ~isempty(additionalArrayArgs)
   if all(cellfun(@isnumeric,additionalArrayArgs)) && ismember(params.inputOutputType,{'3D Array','4D Array'}) %if all arguments are numeric and the input type is Array
      additionalArrayInputs = cellfun(@(x)repmat(x,[size(overlayData{1}) 1]),additionalArrayArgs,'UniformOutput',false); %convert additional arguments to arrays
   else %if any additional argument is not a number
      additionalArrayInputs = cellfun(@(x)num2cell(repmat(x,[size(overlayData{1}) 1])),additionalArrayArgs,'UniformOutput',false); %convert additional arguments to cell arrays
      params.inputOutputType = 'Scalar';  %and force scalar
   end
else
   additionalArrayInputs = {};
end

%parse other additional inputs
additionalArgs = parseArguments(params.additionalArgs,',');

%convert overlays to cell arrays if scalar function
if strcmp(params.inputOutputType,'Scalar')
   for iScan = 1:nScans
      for iOverlay = 1:length(overlayList)
         overlayData{iScan,iOverlay} = num2cell(overlayData{iScan,iOverlay}); %convert overlays to cell arrays
      end
   end
end

%reshape the input if function run separately on each input
if strcmp(params.combinationMode,'Apply function to each overlay')
   overlayData = cellReshape(overlayData,size(overlayData,1),size(overlayData,3),size(overlayData,2));
end


if strcmp(params.combinationMode,'Recursively apply to overlay pairs') %should add additional non-array arguments also ?
   nTotalargs = size(overlayData,2)+length(additionalArrayArgs);
   combineFunctionString = '@(x1';
   for iInput = 2:nTotalargs
      combineFunctionString = [combineFunctionString ',x' num2str(iInput)];
   end
   combineFunctionString = [combineFunctionString ')' ];
   for iInput = nTotalargs:-1:2
      combineFunctionString = [combineFunctionString params.combineFunction '(x' num2str(iInput) ','];
   end
   combineFunctionString = [combineFunctionString 'x1'];
   for iInput = 1:nTotalargs-1
      combineFunctionString = [combineFunctionString ')'];
   end
else
   combineFunctionString = params.combineFunction;
end
combineFunctionHandle = str2func(combineFunctionString);


%-----------------------------construct function call
%output arguments
if params.nOutputOverlays
  functionString = '[';
  for iOutput = 1:params.nOutputOverlays
     functionString = [functionString 'outputData{iScan,' num2str(iOutput) ',iOperations},'];
  end
  functionString(end:end+3)='] = ';
else
  functionString='';
end

%array operator
switch(params.inputOutputType)
  case 'Scalar' %if function operates on scalar
   functionString = [functionString 'cellfun(combineFunctionHandle,']; %we'll use cellfun to apply the function to each voxel
  case '4D Array'                                       %if it operates on 4D arrays
   functionString = [functionString 'feval(combineFunctionHandle,cat(4,']; %we'll concatenate the inputs
  otherwise                                       %if it operates on arrays
   functionString = [functionString 'feval(combineFunctionHandle,'];
end

%input arguments
for iInput = 1:size(overlayData,2)
   functionString = [functionString 'overlayData{iScan,' num2str(iInput) ',iOperations},'];
end
%additional arguments
for iInput = 1:length(additionalArrayInputs)
   functionString = [functionString 'additionalArrayInputs{' num2str(iInput) '},'];
end
if strcmp(params.inputOutputType,'4D Array')
  functionString(end:end+1) = '),'; %add closing bracket to end function cat
end

%additional scalar arguments
for iInput = 1:length(additionalArgs)
   functionString = [functionString 'additionalArgs{' num2str(iInput) '},'];
end

if params.passView
  functionString = [functionString 'thisView,'];
end

functionString(end:end+1) = ');'; %replace the last comma by a closing bracket to end the function



%evaluate the function 
outputData = cell(size(overlayData,1),params.nOutputOverlays, size(overlayData,3));
outputOverlayNames = cell(params.nOutputOverlays, size(overlayData,3));
for iOperations = 1:size(overlayData,3)
  for iScan = 1:size(overlayData,1)
    %check for empty overlays (only if 3D array or scalar)
    emptyInput=false;
    if ismember(params.inputOutputType,{'3D Array','4D Array','Scalar'})
      for iInput = 1:size(overlayData,2)
        if isempty(overlayData{iScan,iInput,iOperations})
          emptyInput=true;
        end
      end
    end
    if emptyInput
      for iOutput = 1:params.nOutputOverlays
        outputData{iScan,iOutput,iOperations} = [];
      end
    else
      try
%         tic
         eval(functionString);
%          toc
      catch exception
        mrWarnDlg(sprintf('There was an error evaluating function %s:\n%s',combineFunctionString,getReport(exception,'basic')));
        set(viewGet(thisView,'figNum'),'Pointer','arrow');drawnow;
        return
      end
      for iOutput = 1:params.nOutputOverlays
    %          if  strcmp(params.inputOutputType,'Scalar')
    %             %convert back to numerical value
    %             outputData{iScan,iOverlay,iOperations} = cell2mat(outputData{iScan,iOverlay,iOperations}); 
    %          end
         %add 0 to all the results to convert logical to doubles, because mrLoadRet doesn't like logical overlays
        switch(params.inputOutputType)
          case 'Structure'
            for jScan = 1:nScans
              outputData{iScan,iOutput,iOperations}.data{jScan} = outputData{iScan,iOutput,iOperations}.data{jScan}+0;
            end
          case {'3D Array','4D Array','Scalar'}
           outputData{iScan,iOutput,iOperations} = outputData{iScan,iOutput,iOperations}+0;
        end
        %check that the size is compatible
        if ~isequal(size(outputData{iScan,iOutput,iOperations}),size(overlayData{iScan}))
          mrWarnDlg(['(combineTransformOverlays) Dimensions of result are not compatible with overlay ([' num2str(size(outputData{iScan,iOutput,iOperations})) '] vs [' num2str(size(overlayData{iScan})) '])']);
        end
      end
    end
  end
end

if params.nOutputOverlays 
  %name of output overlays
  for iOutput=1:params.nOutputOverlays
     if params.nOutputOverlays>1
        name = ['Ouput ' num2str(iOutput) ' - '];
     else
        name = '';
     end
     if ~isempty(params.outputName)
        name = [name params.outputName '('];
     else
        name = [name params.combineFunction '('];
     end
     for iOperations = 1:size(overlayData,3)
        if size(overlayData,3)>1
           outputOverlayNames{iOutput,iOperations} = [name overlayNames{iOperations} ','];
        else
           outputOverlayNames{iOutput,iOperations} = name;
           for iInput = 1:size(overlayData,2)
              outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} overlayNames{iInput} ','];
           end
        end
        for iInput = 1:length(additionalArrayArgs)
           if isnumeric(additionalArrayArgs{iInput})
              outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} num2str(additionalArrayArgs{iInput}) ','];
           else
              outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} additionalArrayArgs{iInput} ','];
           end
        end
        for iInput = 1:length(additionalArgs)
           if isnumeric(additionalArgs{iInput})
              outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} mat2str(additionalArgs{iInput}) ','];
           elseif isa(additionalArgs{iInput},'function_handle')
              outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} func2str(additionalArgs{iInput}) ','];
           else
              outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations} additionalArgs{iInput} ','];
           end
        end
        outputOverlayNames{iOutput,iOperations} = [outputOverlayNames{iOutput,iOperations}(1:end-1) ')'];
     end
  end

  if size(overlayData,3)>1 %reshape the overlay cell array
     outputData = cellReshape(outputData,size(outputData,1),size(outputData,2)*size(outputData,3));
     outputOverlayNames = cellReshape(outputOverlayNames,numel(outputOverlayNames),1);
  end

  %pre-compute coordinates map to put values back from base space to overlay space
  if params.baseSpace && ~params.exportToNewGroup && any(any((base2scan - eye(4))>1e-6)) 
    overlayIndexMap=cell(nScans,1);
    overlayCoordsMap=cell(nScans,1);
    baseCoordsOverlay=cell(nScans,1);
    for iScan=1:nScans
      %make a coordinate map of which overlay voxel each base map voxel corresponds to (convert base coordmap to overlay coord map)
      baseCoordsMap{iScan} = reshape(baseCoordsMap{iScan},numel(baseCoordsMap{iScan})/3,3);
      overlayCoordsMap{iScan} = (base2scan*[baseCoordsMap{iScan}';ones(1,size(baseCoordsMap{iScan},1))])';
      overlayCoordsMap{iScan} = overlayCoordsMap{iScan}(:,1:3);
      overlayCoordsMap{iScan}(all(~overlayCoordsMap{iScan},2),:)=NaN;
      overlayCoordsMap{iScan} = round(overlayCoordsMap{iScan});
      scanDims = viewGet(thisView,'dims',iScan);
      overlayCoordsMap{iScan}(any(overlayCoordsMap{iScan}>repmat(scanDims,size(overlayCoordsMap{iScan},1),1)|overlayCoordsMap{iScan}<1,2),:)=NaN;
      %convert overlay coordinates to overlay indices for manipulation ease
      overlayIndexMap{iScan} = sub2ind(scanDims, overlayCoordsMap{iScan}(:,1), overlayCoordsMap{iScan}(:,2), overlayCoordsMap{iScan}(:,3));

      %now make a coordinate map of which base map voxels each overlay index corresponds to
      %(there will be several maps because each overlay voxels might correspond to several base voxels)

  % %       %METHOD 1
  % %       %sort base indices
  % %       [sortedOverlayIndices,whichBaseIndices] = sort(overlayIndexMap{iScan});
  % %       %remove NaNs (which should be at the end of the vector)
  % %       whichBaseIndices(isnan(sortedOverlayIndices))=[];
  % %       sortedOverlayIndices(isnan(sortedOverlayIndices))=[];
  % %       %find the first instance of each unique index
  % %       firstInstances = sortedIndices(1:end-1) ~= sortedIndices(2:end);
  % %       firstInstances = [true;firstInstances];
  % %       %get the unique overlay indices
  % %       uniqueOverlayIndices = sortedOverlayIndices(firstInstances);
  % %       %compute the number of instances for each  unique overlay index (= number
  % %       %of base different indices for each unique overlay index)
  % %       numberInstances = diff(find([firstInstances;true]));
  % %       maxInstances = max(numberInstances);
  % %       baseCoordsOverlay2{iScan} = sparse(prod(scanDims),maxInstances);
  % %       hWaitBar = mrWaitBar(-inf,'(combineTransformOverlays) Creating base coordinates overlay map for scan');
  % %       %for each unique overlay index, find all the corresponding base indices
  % %       for i = 1:length(uniqueOverlayIndices)
  % %         mrWaitBar( i/length(uniqueOverlayIndices), hWaitBar);
  % %         theseBaseIndices = whichBaseIndices(sortedOverlayIndices==uniqueOverlayIndices(i));
  % %         baseCoordsOverlay2{iScan}(uniqueOverlayIndices(i),1:length(theseBaseIndices))=theseBaseIndices';
  % %       end
  % %       mrCloseDlg(hWaitBar);

      %METHOD 2 (faster)
      %first find the maximum number of base voxels corresponding to a single overlay voxel (this is modified from function 'unique')
      %sort base non-NaN indices
      sortedIndices = sort(overlayIndexMap{iScan}(~isnan(overlayIndexMap{iScan})));
      %find the first instance of each unique index
      firstInstances = sortedIndices(1:end-1) ~= sortedIndices(2:end);
      firstInstances = [true;firstInstances];
      %compute the number of instances for each unique overlay index 
      %(= number of base different indices for each unique overlay index)
      numberInstances = diff(find([firstInstances;true]));
      maxInstances = max(numberInstances);
      baseCoordsOverlay{iScan} = sparse(prod(scanDims),maxInstances);
      %Now for each set of unique overlay indices, find the corresponding base indices
      hWaitBar = mrWaitBar(-inf,'(combineTransformOverlays) Creating base coordinates overlay map for scan');
      for i=1:maxInstances
        mrWaitBar( i/maxInstances, hWaitBar);
        %find set of unique instances of overlay indices
        [uniqueOverlayIndices, whichBaseIndices]= unique(overlayIndexMap{iScan});
        %remove NaNs
        whichBaseIndices(isnan(uniqueOverlayIndices))=[];
        uniqueOverlayIndices(isnan(uniqueOverlayIndices))=[];
        %for each overlay voxel found, set the corresponding base index
        baseCoordsOverlay{iScan}(uniqueOverlayIndices,i)=whichBaseIndices;
        %remove instances that were found from the overlay index map before going through the loop again
        overlayIndexMap{iScan}(whichBaseIndices)=NaN;
      end
      mrCloseDlg(hWaitBar);
    end
  end

  %save the output
  if params.exportToNewGroup
    baseName = [viewGet(thisView,'currentBaseName') 'Volume'];
    groupName=baseName;
    if ismember(groupName,viewGet(thisView,'groupNames'))
      fprintf('(combineTransformOverlays) Installing output overlays in group %s.\n',groupName);
      thisView = viewSet(thisView,'currentGroup',groupName);
      thisView = viewSet(thisView,'currentBase',viewGet(thisView,'baseNum',groupName));
    else
      fprintf('(combineTransformOverlays) Creating group %s.\n',groupName);
      thisView = viewSet(thisView,'newGroup',groupName);
      thisView = viewSet(thisView,'currentGroup',groupName);
      % export the flat map as an empty volume
      base = viewGet(thisView,'baseCache');
      baseNum = viewGet(thisView,'currentBase');
      if isempty(base)
        base.im = getBaseSlice(thisView,viewGet(thisView,'curslice'),viewGet(thisView,'baseSliceIndex',baseNum),viewGet(thisView,'rotate'),baseNum,viewGet(thisView,'basetype'));
      end
      baseVolume = viewGet(thisView,'baseVolume');
      hdr = baseVolume.hdr;
      hdr.bitpix = 32;   
      hdr.datatype = 16;
      hdr.is_analyze = 1;
      hdr.scl_slope = 1;
      hdr.endian = 'l';
      voxelSize(1:2) = repmat(mean(voxelSize(1:2)),1,2);
      if any(voxelSize ~= viewGet(thisView,'basevoxelsize',baseNum))
         hdr.pixdim = [0 voxelSize 0 0 0 0]';        % all pix dims must be specified here
         hdr.qform44 = diag([voxelSize 0]);
         hdr.sform44 = hdr.qform44;
      end
      % Copy file to the tseries directory
      tseriesDir = viewGet(thisView,'tseriesDir');
      scanFileName = [baseName mrGetPref('niftiFileExtension')];
      newPathStr = fullfile(tseriesDir,scanFileName);
      [bytes,hdr] = cbiWriteNifti(newPathStr,repmat(base.im,[1 1 size(outputData{1},3)]),hdr);
      % Add it
      scanParams.fileName = scanFileName;
      thisView = viewSet(thisView,'newScan',scanParams);
      % create dummy analysis
      analysis.name = 'combineTransformOverlays';
      fprintf('(combineTransformOverlays) Creating Analysis %s.\n',analysis.name);
      analysis.type = 'dummy';
      analysis.groupName = baseName;
      analysis.function = '';
      analysis.reconcileFunction = 'dummyAnalysisReconcileParams';
      analysis.guiFunction = '';
      analysis.params = [];
      analysis.overlays =[];
      analysis.curOverlay = [];
      analysis.date = datestr(now);
%       analysis.clipAcrossOVerlays = 0; %this doesn't work
      thisView = viewSet(thisView,'newanalysis',analysis);
      thisView = viewSet(thisView,'clipAcrossOverlays',0); %this doesn't work
      
      %if the base is a flat map, need to create a new volume flat map
      if viewGet(thisView,'basetype')==1
        thisView = loadAnat(thisView,getLastDir(newPathStr),fileparts(newPathStr));
        saveAnat(thisView,getLastDir(newPathStr));
        thisView.baseVolumes(thisView.curBase).clip=[0 1]; %should add a viewSet case for this but I don't have time
        thisView = viewSet(thisView,'rotate',0); %this doesn't work
      end
    end
  else
    groupName=viewGet(thisView,'groupName');
  end
  
  defaultOverlay.groupName = groupName;
  defaultOverlay.function = 'combineTransformOverlays';
  defaultOverlay.interrogator = '';
  defaultOverlay.type = 'combination';
  defaultOverlay.params = params;
  defaultOverlay.clip = [];
  defaultOverlay.range = [];
  defaultOverlay.name = [];
  defaultOverlay.data = [];
  for iOverlay = 1:size(outputData,2)
    switch(params.inputOutputType)
      case {'3D Array','4D Array','Scalar'}
        outputOverlay(iOverlay) = defaultOverlay;
        outputOverlay(iOverlay).data = outputData(:,iOverlay);
        for iOutput = 1:size(outputData,1)
          isNotEmpty(iOutput) = ~isempty(outputData{iOutput,iOverlay});
        end
        allScansData = cell2mat(outputData(isNotEmpty,iOverlay));
      case 'Structure'
        outputOverlay(iOverlay) = copyFields(defaultOverlay,outputData{iOverlay});
        allScansData = cell2mat(outputOverlay(iOverlay).data);
    end
    if ~params.exportToNewGroup && params.baseSpace && any(any((base2scan - eye(4))>1e-6)) %put back into scan/overlay space
      for iScan=1:nScans
        if ~isempty(outputOverlay(iOverlay).data{iScan}) 
          if viewGet(thisView,'basetype')==1
            data = zeros(scanDims);
            datapoints=zeros(prod(scanDims),1);
            for i=1:size(baseCoordsOverlay{iScan},2)
              thisBaseCoordsMap = full(baseCoordsOverlay{iScan}(:,i));
              data(logical(thisBaseCoordsMap)) = data(logical(thisBaseCoordsMap)) + ...
                    outputOverlay(iOverlay).data{iScan}(thisBaseCoordsMap(logical(thisBaseCoordsMap)));
              datapoints = datapoints+logical(thisBaseCoordsMap);
            end
            datapoints = reshape(datapoints,scanDims);
            outputOverlay(iOverlay).data{iScan} = data ./datapoints;
          else
            keyboard %not implemented
          end
        end
      end
    end
    maxValue = max(allScansData(allScansData<inf));
    minValue = min(allScansData(allScansData>-inf));
    outputOverlay(iOverlay).clip = [minValue maxValue];
    outputOverlay(iOverlay).range = [minValue maxValue];
    outputOverlay(iOverlay).name = outputOverlayNames{iOverlay};
  end
  thisView = viewSet(thisView,'newoverlay',outputOverlay);

  refreshMLRDisplay(thisView.viewNum);
end
set(viewGet(thisView,'figNum'),'Pointer','arrow');drawnow;

function [arguments, nArgs] = parseArguments(argumentString, separator)
   
%parse string of arguments separated by separator and put them into a cell array of numerical and string arguments
%non-numerical values that are not between quotes are converted into strings
%
% Julien Besle, 08/07/2010
nArgs = 0;
arguments = cell(0);
remain = argumentString;
while ~isempty(remain)
   nArgs = nArgs+1;
   [token,remain] = strtok(remain, separator);
   try
      arguments{nArgs} = eval(token);
   catch exception
      if strcmp(exception.identifier,'MATLAB:UndefinedFunction')
         arguments{nArgs} = token;
      else
         mrErrorDlg(['(parseArguments) could not read argument: ' exception.message]);
      end
   end
      
end

function printHelp(params)

if strcmp(params.combineFunction,'User Defined')
  mrWarnDlg('(combineTransformOverlays) Please select a combination function');
else
  helpString = help(params.combineFunction);
  if isempty(helpString)
    mrWarnDlg('(combineTransformOverlays) No help available');
  else
    fprintf('\n');
    disp(helpString);
  end
end

