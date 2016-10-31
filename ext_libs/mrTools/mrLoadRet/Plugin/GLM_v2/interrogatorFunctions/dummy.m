% dummyInterrogator 
%
%        $Id$
%      usage: [  ] = dummyInterrogator(thisView,overlayNum,scanNum,x,y,z,roi)
%         by: julien besle
%       date: 2010-03-09
%     inputs: 
%    outputs: 
%
%    purpose: easy access to currently loaded data

function dummyInterrogator(thisView,overlayNum,scanNum,x,y,z,roi)

keyboard

%get the current overlay data
overlay = viewGet(thisView,'overlay',overlayNum);

%number of values < .05
disp([ thisView.analyses{thisView.curAnalysis}.name ': ' overlay.name ': values less than .05: ' num2str(length(find(overlay.data{1}<.05)) / length(find(~isnan(overlay.data{1})))) ]);


return;
