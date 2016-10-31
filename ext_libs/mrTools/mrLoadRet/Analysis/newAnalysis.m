% newAnalysis.m
%
%        $Id: newAnalysis.m 1942 2010-12-16 18:14:41Z julien $ 
%      usage: v = newAnalysis(v,analysisName)
%         by: justin gardner, taken out from mrLoadRetGUI by julien besle
%       date: 2014/06/06
%    purpose: creates a new empty analysis

function v = newAnalysis(v,analysisName)

analysis.name = analysisName;
analysis.type = 'dummy';
analysis.groupName = viewGet(v,'groupName',viewGet(v,'currentGroup'));
analysis.function = 'dummyAnalysis';
analysis.reconcileFunction = 'dummyAnalysisReconcileParams';
analysis.reconcileFunction = 'dummyAnalysisMergeParams';
analysis.guiFunction = 'dummyAnalysisGUI';
analysis.params = [];
analysis.overlays =[];
analysis.curOverlay = [];
analysis.date = datestr(now);
v = viewSet(v,'newanalysis',analysis);
