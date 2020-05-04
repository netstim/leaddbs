function ea_ttest(x,y,description,ylabel,groups,colors)
% this simple function is a small wrapper for a t-test figure and has been
% written to go in line with the ea_corrplot wrapper.

if ~exist('groups', 'var')
    groups = {'0', '1'};
end

if exist('colors', 'var') && ~isempty(colors)
    map = colors;
    if isnumeric(map)
        if size(map,1) ~= 2
            error('Number of colors should be 2 for two sample t-test!');
        end
    end
else
    map = [0.4127 0.4127 1;1 0.4127 0.4127];
end

if ischar(map)
    colorOptions = {'map', map};
else
    colorOptions = {'map', map, 'n_color', size(map,1), 'n_lightness', 1};
end

disp(description);

[h,p,ci,stats] = ttest2(x,y)

data = [x;y]; % concat for gramm
colorVar = groups([ones(length(x),1);repmat(2,length(y),1)]);
g(1) = gramm('x',ones(length(data),1),'y',data,'color',colorVar);
g(1).set_color_options(colorOptions{:});
g(1).stat_boxplot();
g(1).set_title(description);

g.set_names('x','','y',ylabel,'color','Group');
g(1).axe_property('XTickLabel','');
g(1).axe_property('XTick',[]);

figure('Position',[100 100 360 600]);
g.draw();
