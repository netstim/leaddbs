function params = motionCompGUI(varargin);
% params = averageTSeries('groupName',groupName,'params',params)
% params = averageTSeriesGUI('groupName',groupName);
% params = averageTSeriesGUI('params',params);
% 
% This will create a GUI for motion compensation, but for now it is just a
% hack specifying default parameters.

% Parse varargin
for index = 1:2:length(varargin)
    field = varargin{index};
    val = varargin{index+1};
    switch field
        case 'groupName'
            groupName = val;
        case 'params'
            params = val;
        otherwise
            mrWarnDlg('Invalid initialization argument')
    end
end

% Reconcile/initialize params with current status of group and ensure that
% it has the required fields.
if ieNotDefined('params')
    params = motionCompReconcileParams(groupName);
else
    params = motionCompReconcileParams(groupName,params);
end
