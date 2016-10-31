function [value prefDefaults] = mrGetPref(pref)
%
% value = mrGetPref(pref)
%
% Replaces Matlab's getpref function. Gets a field from the global variable
% mrDEFAULTS.preferences, which is a structure with fields for each
% preference. Returns the value for that preference.
%
% Examples:
%   mrGetPref('verbose');
%   mrGetPref('site');
%   mrGetPref('niftiFileExtension');
%   mrGetPref('interpMethod');
%   mrGetPref('overwritePolicy');
%   mrGetPref('maxBlockSize');
%   mrGetPref('volumeDirectory');
%
% Note that the mrDefaults file is usually saved in ~/.mrDefaults
% but that location can be overridden (see mrDefaultsFilename.m)
%
%   mrGetPref with no arguments returns a list of known preference and known default values
%   [prefNames prefDefaults] = mrGetPref;
%
%   if you type mrGetPref alone, it will print out all known preferences
%   mrGetPref;
%
% djh, 5/2007
% %	$Id: mrGetPref.m 2890 2013-10-27 12:52:21Z justin $	

% with no arguments, return a list of possible preferences
prefNames = {'overwritePolicy','verbose','graphWindow','checkParamsConsistency'...
   'maxBlocksize','roiCacheSize','baseCacheSize','overlayCacheSize','defaultPrecision',...
   'defaultInterrogators','systemInterrogators',...,
   'importROIPath','volumeDirectory','niftiFileExtension','fslPath',...
   'selectedROIColor','roiContourWidth','roiPolygonMethod',...,
   'interpMethod','corticalDepthBins','roiCorticalDepthDisplayRatio', 'multiSliceProjectionMethod','colorBlending',...
   'pluginPaths','selectedPlugins',...
   'statisticalTestOutput',...
   'site','magnet','coil','pulseSequence',...
   'maxArrayWidthForParamsDialog','maxArrayHeightForParamsDialog',...
   'mlrVolDisplayControls','mlrVolOverlayAlpha','motionCompDefaultParams','colorNames',...
    'mlrPath','vistaPath','lastPath'...
	    };

% set the defaults for preference we have defaults for. Note that the "find" in
% here is to make sure that the prefDefaults list matches the prefNames order
prefDefaults{length(prefNames)} = [];
prefDefaults{find(strcmp('overwritePolicy',prefNames))} = {'Ask','Merge','Rename','Overwrite'};
prefDefaults{find(strcmp('verbose',prefNames))} = {'No','Yes'};
prefDefaults{find(strcmp('graphWindow',prefNames))} = {'Replace','Make new'};
prefDefaults{find(strcmp('checkParamsConsistency',prefNames))} = {'Yes','No'};
prefDefaults{find(strcmp('maxBlocksize',prefNames))} = 250000000;
prefDefaults{find(strcmp('roiCacheSize',prefNames))} = 100;
prefDefaults{find(strcmp('baseCacheSize',prefNames))} = 50;
prefDefaults{find(strcmp('overlayCacheSize',prefNames))} = 50;
prefDefaults{find(strcmp('defaultPrecision',prefNames))} = 'double';
prefDefaults{find(strcmp('volumeDirectory',prefNames))} = '';
prefDefaults{find(strcmp('niftiFileExtension',prefNames))} = {'.img','.nii'};
prefDefaults{find(strcmp('fslPath',prefNames))} = 'FSL not installed';
prefDefaults{find(strcmp('selectedROIColor',prefNames))} = color2RGB;
prefDefaults{find(strcmp('selectedROIColor',prefNames))}{end+1} = 'none';
prefDefaults{find(strcmp('roiContourWidth',prefNames))} = 1;
prefDefaults{find(strcmp('roiCorticalDepthDisplayRatio',prefNames))} = .5;
prefDefaults{find(strcmp('roiPolygonMethod',prefNames))} = {'getpts','roipoly','getptsNoDoubleClick'};
prefDefaults{find(strcmp('interpMethod',prefNames))} = {'nearest','linear','spline','cubic'};
prefDefaults{find(strcmp('corticalDepthBins',prefNames))} = 11;
prefDefaults{find(strcmp('multiSliceProjectionMethod',prefNames))} = {'Average','Maximum Intensity Projection'};
prefDefaults{find(strcmp('colorBlending',prefNames))} = {'Additive','Alpha blend'};
prefDefaults{find(strcmp('pluginPaths',prefNames))} = '';
prefDefaults{find(strcmp('selectedPlugins',prefNames))} = '';
prefDefaults{find(strcmp('statisticalTestOutput',prefNames))} = {'P value','Z value','-log10(P) value'};
prefDefaults{find(strcmp('site',prefNames))} = 'NYU';
prefDefaults{find(strcmp('magnet',prefNames))} = {{'Allegra 3T','other'}};
prefDefaults{find(strcmp('coil',prefNames))} = {{'LifeService','Siemens birdcage','Nova birdcage','Nova surface','Nova quadrapus','Nova visual array','other'}};
prefDefaults{find(strcmp('pulseSequence',prefNames))} = {{'cbi_ep2d_bold','other'}};
prefDefaults{find(strcmp('maxArrayWidthForParamsDialog',prefNames))} = 25;
prefDefaults{find(strcmp('maxArrayHeightForParamsDialog',prefNames))} = 50;
prefDefaults{find(strcmp('mlrVolDisplayControls',prefNames))} = false;
prefDefaults{find(strcmp('mlrVolOverlayAlpha',prefNames))} = 0.8;
prefDefaults{find(strcmp('motionCompDefaultParams',prefNames))} = [];
prefDefaults{find(strcmp('colorNames',prefNames))} = {};
prefDefaults{find(strcmp('mlrPath',prefNames))} = '';
prefDefaults{find(strcmp('vistaPath',prefNames))} = '';
prefDefaults{find(strcmp('lastPath',prefNames))} = '';

if nargin == 0
  if nargout > 0
    % return arguments
    value = prefNames;
  else
    % print out list of preferences
    for i = 1:length(prefNames)
      %get the preference
      prefValue =  mrGetPref(prefNames{i});
      % print it out
      if isnumeric(prefValue)
	disp(sprintf('%s: %s',prefNames{i},num2str(prefValue)));
      elseif isstr(prefValue)
	disp(sprintf('%s: %s',prefNames{i},prefValue));
      elseif iscell(prefValue)
	mrDisp(sprintf('%s:',prefNames{i}));
	for j = 1:length(prefValue)
	  if isnumeric(prefValue{j})
	    mrDisp(sprintf(' %s',num2str(prefValue{j})));
	  elseif isstr(prefValue{j})
	    mrDisp(sprintf(' ''%s''',prefValue{j}));
	  end
	end
	mrDisp(sprintf('\n'));
      end
    end
  end
  return
end

global mrDEFAULTS

% fix any caps differences
prefNum = find(strcmp(lower(pref),lower(prefNames)));
if ~isempty(prefNum)
  pref = prefNames{prefNum};
end

% read the preferences and figlocs
if isempty(mrDEFAULTS)
  mrDEFAULTS = loadMrDefaults;
end

if isfield(mrDEFAULTS.prefs,pref)
    value = getfield(mrDEFAULTS.prefs,pref);
else
  % not set yet, take the top most possibility in the default
  % list, otherwise return empty
  if ~isempty(prefNum) && ~isempty(prefDefaults{prefNum})
    if iscell(prefDefaults{prefNum})
      value = prefDefaults{prefNum}{1};
    else
      value = prefDefaults{prefNum};
    end
  else
    value = [];
  end
end

% default value for selectedROIColor
if strcmp(pref,'selectedROIColor') && isempty(value)
  value = 'white';
end

% deal with interrogators
if strcmp(pref,'systemInterrogators')
  if isempty(value)
    value = {'timecoursePlot','makeFlat','searchForVoxel'};
  end
end
if strcmp(pref,'defaultInterrogators')
  if isempty(value)
    value = mrGetPref('systemInterrogators');
  else
    value = union(value,mrGetPref('systemInterrogators'));
  end
end
