% eventRelatedPreProcess.m
%
%        $Id$
%      usage: d = eventRelatedPreProcess(d)
%         by: justin gardner
%       date: 03/29/07
%    purpose: make various changes to the d structure before
%             passing it on to r2 processing, this allows for
%             little modifications that might be necessary to
%             fix small things
%
function [d, failure] = eventRelatedPreProcess(d,type,verbose)

% check arguments
if ~ismember(nargin,[2 3])
  help eventRelatedPreProcess
  return
end
failure = 0;
if ieNotDefined('verbose')
   verbose = 1;
end

% nothing to do
if ieNotDefined('type'),return,end

% keep track of whether this is a function name
isRunFunctionName = 1;

% move all the stimvolumes by some amount
stimvolShift  = getfieldnum(type,'stimvolShift');
if ~isempty(stimvolShift)
  disp(sprintf('Shifting by %i volumes',stimvolShift));
  for stimnum = 1:length(d.stimvol)
    d.stimvol{stimnum} = d.stimvol{stimnum}+stimvolShift;
  end
  isRunFunctionName = 0;
end

% length of tr
tr = getfieldnum(type,'tr');
if ~isempty(tr)
  disp(sprintf('Setting TR to %0.2f',tr));
  d.tr = tr;
  isRunFunctionName = 0;
end

% preprocessing for sdt experiment
runName = getfieldstr(type,'run');
% or look for a type that is just the name of a function
if ~isempty(type) && isstr(type) && exist(type,'file')
  runName = type;
end

if isRunFunctionName && isempty(runName)
  mrWarnDlg(sprintf('(eventRelatedPreProcess) Could not run pre process function %s',type));
  failure = 1;
end

% if we have a runName then run the preprocess
if ~isempty(runName)
  % get run functiona and arguments
  [runFun runName] = strtok(runName,'()');
  [runArgs] = strtok(runName,'()');
  % default to passing the d
  if isempty(runArgs),runArgs = 'd';end
  % now see if it is a function and if so, run it
  if exist(sprintf('%s.m',stripext(runFun)))==2
     if verbose
         disp(sprintf('Running d=%s(%s)',runFun,runArgs));
     end
    d = eval(sprintf('%s(%s)',runFun,runArgs));
  else
    disp(sprintf('(eventRelatedPreProcess) Could not find preprocess script %s',runFun));
    failure = 1;
    keyboard
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if the type has something like fieldname=num
% it returns the num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fieldval = getfieldnum(type,fieldname)

fieldval = getfieldstr(type,fieldname);
fieldval = str2num(fieldval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if the type has something like fieldname=string
% it returns the string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fieldval = getfieldstr(type,fieldname)

% see if it is there
if ~strfind(type,sprintf('%s=',fieldname))
  fieldval = [];
end

% get the field
fieldval = type(first(strfind(type,sprintf('%s=',fieldname))):length(type));
% up to the end
fieldval = fieldval(first(strfind(fieldval,'='))+1:length(fieldval));
% remove everything after first space
if (strfind(fieldval,' '))
  fieldval = fieldval(1:first(strfind(fieldval,' '))-1);
end


