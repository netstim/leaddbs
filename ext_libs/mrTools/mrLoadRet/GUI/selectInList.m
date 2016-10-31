% selectInList.m
%       
%      usage: [numberList,names] = selectInList(thisView,type,title,preselected,groupNum)
%         by: julien besle
%       date: 20/12/2010
%        $Id$
%    purpose: asks users to select a list of a given type of mrLoadRet variables from the thisView
%                 type can be: 'overlays','scans','analyses','bases','rois', 'groups



function [numberList,names] = selectInList(thisView,type,title,preselected,groupNum)

switch(lower(type))
  case {'scans','scan'}
    type='Scans';
    if ieNotDefined('groupNum')
       groupNum = viewGet(thisView,'currentGroup');
    end
    numberInView = viewGet(thisView,'nScans',groupNum);
    for i = 1:numberInView
      names{i} = sprintf('%i:%s (%s)',i,viewGet(thisView,'description',i,groupNum),viewGet(thisView,'tSeriesFile',i,groupNum));
    end
    
  case {'overlays','overlay'}
    type='Overlays';
    names = viewGet(thisView,'overlayNames');
    
  case {'analyses','analysis'}
    type='Analysis';
    names = viewGet(thisView,'analysisNames');
    
  case {'bases','base'}
    type='Bases';
    names = viewGet(thisView,'baseNames');
    
  case {'rois','roi'}
    type='ROIs';
    names = viewGet(thisView,'roiNames');
    
  case {'groups','group'}
    type='Groups';
    names = viewGet(thisView,'groupNames');
    
  otherwise
    mrWarnDlg(['(selectInList) Unknown type ''' type '''']);
    numberList=[];
    names=[];
    return;

end
numberInView = length(names);

%Check for zero:
if numberInView == 0
  mrWarnDlg(['(selectInList) No ' type ' found!']);
  numberList=[];
  names=[];
  return
end


if ieNotDefined('title')
  title = ['Choose ' type];
end
if ~exist('preselected','var') %if preselected exists but is empty, leave as it is
  switch(lower(type))
    case {'scans','scan'}
      preselected = viewGet(thisView,'currentScan');

    case {'overlays','overlay'}
      preselected = viewGet(thisView,'currentOverlay');

    case {'analyses','analysis'}
      preselected = viewGet(thisView,'currentAnalysis');

    case {'bases','base'}
      preselected = viewGet(thisView,'currentBase');

    case {'rois','roi'}
      preselected = viewGet(thisView,'currentROI');
      
    case {'groups','group'}
      preselected = viewGet(thisView,'currentGroup');
  end
end

%add a line for all
if numberInView>1
  names = [names {'All'}];
  numberInView = numberInView+1;
end

preselection = zeros(1,numberInView);
if ~nnz(~preselected)
  preselection(preselected) = 1;
end

iSel = buttondlg(title, names,preselection);
if isempty(iSel)
  numberList = iSel; %if cancel has been pressed, this will be a 0*0 matrix, 
  %but if the top close button has been pressed, it will be a 0*1 matrix
else
  if numberInView>1 && iSel(end) %if 'All' has been selected
    numberList = 1:numberInView-1;
  else
    numberList = find(iSel); %if OK is pressed but nothing has been selected, this will output a 1*0 array
  end
end

return;
