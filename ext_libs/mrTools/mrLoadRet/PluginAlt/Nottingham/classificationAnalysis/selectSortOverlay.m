function overlay = selectSortOverlay(view)
   groupName = viewGet(view,'groupName');
paramsInfo = {...
    {'groupName',putOnTopOfList(viewGet(view,'groupName'),viewGet(view,'groupNames')),'Name of group from which to make concatenation'}};


x = mrParamsDialog(paramsInfo);

thisView = viewSet(view,'groupName',x.groupName);
analysisNames = viewGet(thisView,'analysisNames');

paramsInfo ={...
    {'analysis',analysisNames,'xxxxx'}};

x = mrParamsDialog(paramsInfo);

a_num = viewGet(thisView,'analysisnum',x.analysis);

overlayNames = viewGet(thisView,'overlayNames',a_num);

paramsInfo ={...
    {'overlays',overlayNames,'xxxxx'}};
x = mrParamsDialog(paramsInfo);

o_num = viewGet(thisView,'overlayNum',x.overlays,a_num);
overlay = viewGet(thisView,'overlay',o_num,a_num);

scan = viewGet(thisView,'nScans');
if scan>1
    paramsInfo={...
        {'scan',scan,'xxx'}};
    x = mrParamsDialog(paramsInfo);
    scan = x.scan;
end

overlay.data = overlay.data(scan);

view = viewSet(view,'groupName',groupName);