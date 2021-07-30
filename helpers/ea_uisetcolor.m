function color = ea_uisetcolor(varargin)
% Wrapper for MATLAB built-in function 'uisetcolor'
%
% The original color picker uses 'jxbrowser-chromium' as backend since
% version R2016b (9.1), which has a bug on Linux. All the MATLAB component
% uses 'jxbrowser-chromium' (like color picker and help browser) might hang
% when closing the window. The color picker is also slower now (on other
% platforms, too). This warpper will set the color picker to the old
% fashion for MATLAB R2016b and later versions on Linux. It can be used
% until MATLAB fixes the bug in future release.
%
% Set to old version for MATLAB R2016b, R2017a, R2017b, R2018a and R2019b+

if ismac || ispc
    color = uisetcolor(varargin{:});
    return;
end

if isMatlabVer('<',[9,1]) ... % < R2016b
        || isMatlabVer('==',[9,5]) ... % R2018b
        || isMatlabVer('==',[9,6]) % R2019a
    color = uisetcolor(varargin{:});
else
    if isMatlabVer('==',[9,1])	% R2016b
        currVer = getpref('Mathworks_uisetcolor', 'Version');
        setpref('Mathworks_uisetcolor', 'Version', 1);
        color = uisetcolor(varargin{:});
        setpref('Mathworks_uisetcolor', 'Version', currVer);
    else	% R2017a, R2017b, R2018a, R2019b+
        newColorPicker = 'matlab.ui.internal.dialog.WebColorChooser';
        oldColorPicker = 'matlab.ui.internal.dialog.ColorChooser';
        s = settings;
        % 'TemporaryValue' only effective in the current MATLAB session
        s.matlab.ui.dialog.uisetcolor.ControllerName.TemporaryValue = oldColorPicker;
        color = uisetcolor(varargin{:});
        s.matlab.ui.dialog.uisetcolor.ControllerName.TemporaryValue = newColorPicker;
    end
end
