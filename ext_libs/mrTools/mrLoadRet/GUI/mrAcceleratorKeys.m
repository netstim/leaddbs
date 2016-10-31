% mrAcceleratorKeys.m
%
%      usage: mrAccerleatorKeys(event,viewNum,val)
%         by: justin gardner
%       date: 06/01/07
%    purpose: set up accelerator keys
%             start by calling
%             mrAcceleratorKeys('init',viewNum);
%             turn off
%             mrAcceleratorKeys('end',viewNum);
function retval = mrAcceleratorKeys(event,viewNum,val)

% check arguments
if ~any(nargin == [1 2 3])
  help mrAcceleratorKeys
  return
end

switch (event)
 case 'init'
  initHandler(viewNum);
 case 'end'
  endHandler(viewNum);
 case 'keyPress'
  keyPressHandler(viewNum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mouseup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keyPressHandler(viewNum)

mrGlobals;
v = MLR.views{viewNum};
currentChar = get(MLR.views{viewNum}.figure,'CurrentCharacter');

if ~isempty(currentChar)
  switch double(currentChar)
   % left arrow
   case 28
    curSlice = viewGet(v,'curSlice');
    viewSet(viewNum,'curslice',curSlice-1);
    refreshMLRDisplay(viewNum);
   % right arrow
   case 29
    curSlice = viewGet(v,'curSlice');
    viewSet(viewNum,'curslice',curSlice+1);
    refreshMLRDisplay(viewNum);
   % up arrow
   case 30
    curScan = viewGet(v,'curScan');
    viewSet(viewNum,'curscan',curScan+1);
    refreshMLRDisplay(viewNum);
   % down arrow
   case 31
    curScan = viewGet(MLR.views{viewNum},'curScan');
    viewSet(viewNum,'curscan',max(1,curScan-1));
    refreshMLRDisplay(viewNum);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end the mrAcceleratorKeys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function endHandler(viewNum)

mrGlobals;
set( MLR.views{viewNum}.figure,'KeyPressFcn','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the interrogator handler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initHandler(viewNum)

mrGlobals;

% set the callbacks appropriately
set( MLR.views{viewNum}.figure,'KeyPressFcn',sprintf('mrAcceleratorKeys(''keyPress'',%i)',viewNum));
