function cg = ea_colorgradient(steps, color1, color2, color3)
% Generate linear color gradient (used for customized colormap)

if ~exist('steps', 'var')   % Input steps
    steps = inputdlg('Please enter the number of the steps:','Steps',1,{'256'});
    steps = str2double(steps);
end

if ~exist('color1', 'var')  % Input color1
    color1 = ea_uisetcolor('Select color 1');
end

if ~exist('color2', 'var')  % Input color2
    color2 = ea_uisetcolor('Select color 2');
end

if nargin <= 2   % Only color1 specified
    color3 = ea_uisetcolor('Select color 3');
    if isscalar(color3) % Cancelled input, color3 == 0
        clear color3
    end
end

if ischar(color1)   % Convert HEX color to RGB values
    color1 = ea_hex2rgb(color1);
end

if ischar(color2)   % Convert HEX color to RGB values
    color2 = ea_hex2rgb(color2);
end

% Generate color gradient
if ~exist('color3', 'var')
    cg = [linspace(color1(1), color2(1), steps)', ...
          linspace(color1(2), color2(2), steps)', ...
          linspace(color1(3), color2(3), steps)'];
else
    if ischar(color3)	% Convert HEX color to RGB values
        color3 = ea_hex2rgb(color3);
    end

    cg1 = [linspace(color1(1), color2(1), round(steps/2))', ...
           linspace(color1(2), color2(2), round(steps/2))', ...
           linspace(color1(3), color2(3), round(steps/2))'];
    cg2 = [linspace(color2(1), color3(1), round(steps/2))', ...
           linspace(color2(2), color3(2), round(steps/2))', ...
           linspace(color2(3), color3(3), round(steps/2))'];
    
    if mod(steps, 2) == 0
        cg = [cg1; cg2];
    else
        cg = [cg1; cg2(2:end,:)];
    end
end
