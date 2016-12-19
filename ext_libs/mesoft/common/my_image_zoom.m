function erg= my_image_zoom(small, high, axesHd, viewStr)
%function erg= my_image_zoom(small, high, axesHd, viewStr)
%
%  Zum richtigem zoomen, so dass bilder und Linien richtig zu einander stehen.
%
%  my_image_zoom([30 22 10], [110 70 10], gca, 'z-axis')
%
% Bjoern W. Kreher
% 02/04
%
% UNIX

erg= [];

if ~exist('axesHd') || isempty(axesHd)
    axesHd= gca;
end

if ~exist('viewStr') || isempty(viewStr)
    viewStr= 'z-axis';
end


childHd= get(axesHd, 'Children');

imagHd= [];
for i= 1:length(childHd)
    if strcmp(get(childHd(i), 'type'), 'image')
        imagHd= childHd(i);
    end
end

if isempty(imagHd)
    return;
end

sizeAy= size(get(imagHd, 'CData'));

erg= sizeAy;

axes(axesHd);
if strcmp(viewStr, 'z-axis')
    axis([small(2) - 0.5 high(2) + 0.5 sizeAy(1) - high(1) + 1 - 0.5 sizeAy(1) - small(1) + 1 + 0.5]);
%    set(gca, 'XLim', [small(2) - 0.5, high(2) + 0.5], 'YLim', sizeAy(1) - [high(1) + 0.5, small(1) - 0.5] + 1);
elseif strcmp(viewStr, 'y-axis')
    axis([small(1) - 0.5 high(1) + 0.5 small(3) - 0.5 high(3) + 0.5])
%    set(gca, 'YLim', [small(3) - 0.5, high(3) + 0.5], 'XLim', [small(1) - 0.5, high(1) + 0.5]);
elseif strcmp(viewStr, 'x-axis')
    axis([small(2) - 0.5 high(2) + 0.5 small(3) - 0.5 high(3) + 0.5])
%    set(gca, 'YLim', [small(3) - 0.5, high(3) + 0.5], 'XLim', [small(2) - 0.5, high(2) + 0.5]);
end
        
        