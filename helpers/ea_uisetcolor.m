function color = ea_uisetcolor(varargin)
% Wrapper for MATLAB built-in function 'uisetcolor'
%
% The original color picker uses 'jxbrowser-chromium' as backend since
% version R2016b (9.1), which has a bug on Linux. All the MATLAB component uses
% 'jxbrowser-chromium' (like color picker and help browser) might hang when
% closing the window. The color picker is also slower now (on other
% platforms, too). This warpper will set the color picker to the old
% fashion for MATLAB R2016b and later versions on Linux. It can be used
% until MATLAB fixes the bug in future release.
%
% Seems fixed on R2018b (9.5). Reocurring in R2019b (9.7)...

ver = version;
ver = str2double(ver(1:3));

if ismac || ispc || ver < 9.1 %|| ver == 9.5
    color = uisetcolor(varargin{:});
else
    if ver == 9.1   % R2016b
        currVer = getpref('Mathworks_uisetcolor', 'Version');
        setpref('Mathworks_uisetcolor', 'Version', 1);
        color = uisetcolor(varargin{:});
        setpref('Mathworks_uisetcolor', 'Version', currVer);
    else    % Since R2017a
        s = settings;
        oldColorPicker = 'matlab.ui.internal.dialog.ColorChooser';
        % For reference:
        % newColorPicker = 'matlab.ui.internal.dialog.WebColorChooser';

        % 'TemporaryValue' only effective in the current MATLAB session
        s.matlab.ui.dialog.uisetcolor.ControllerName.TemporaryValue = oldColorPicker;
        color = uisetcolor(varargin{:});
    end
end
