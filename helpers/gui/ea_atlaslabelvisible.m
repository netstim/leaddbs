function ea_atlaslabelvisible(hobj,~,labels,onoff)
% Match each label to its scene toggle, including missing/non-text labels.
% Atlas presets index structures, whereas labels and colorbuttons index the
% rendered hemispheres; toolbar order is not a reliable mapping either.
labels=labels(:);
isLabel=arrayfun(@(label) isa(label,'matlab.graphics.primitive.Text'),labels);
if isempty(hobj) || strcmp(onoff,'off')
    if any(isLabel)
        set(labels(isLabel),'Visible',onoff);
    end
    return
end

resultfig=ancestor(hobj,'figure');
colorbuttons=getappdata(resultfig,'colorbuttons');
for sceneIndex=find(isLabel)'
    visible='off';
    if sceneIndex<=numel(colorbuttons) && isgraphics(colorbuttons(sceneIndex)) && ...
            strcmp(colorbuttons(sceneIndex).State,'on')
        visible='on';
    end
    labels(sceneIndex).Visible=visible;
end
end
