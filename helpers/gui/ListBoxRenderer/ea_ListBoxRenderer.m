function ea_ListBoxRenderer(control)
% Add dynamic tooltip to ListBox and Popupmenu by overriding the Renderer

% Skip on Windows, since the width of the popup menu will be extended
% automatically if entries inside are too long.
if ispc
    return;
end

if exist('ListBoxRenderer','class') ~= 8
    ea_checkJavaClassPath;
    return;
end

switch control.Style
    case 'listbox'
        jScrollPane = findjobj_fast(control);
        if ~isempty(jScrollPane)
            jScrollPane.getViewport.getView.setCellRenderer(ListBoxRenderer);
        end
    case 'popupmenu'
        jComboBox = findjobj_fast(control);
        if ~isempty(jComboBox)
            jComboBox.setRenderer(ListBoxRenderer);
        end
end
