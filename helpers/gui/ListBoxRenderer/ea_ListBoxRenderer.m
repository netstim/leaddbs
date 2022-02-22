function ea_ListBoxRenderer(control)
% Add dynamic tooltip to ListBox and Popupmenu by overriding Renderer

if exist('ListBoxRenderer','class') ~= 8
    ea_checkJavaClassPath;
    return;
end

switch control.Style
    case 'listbox'
        jScrollPane = findjobj_fast(control);
        jScrollPane.getViewport.getView.setCellRenderer(ListBoxRenderer);
    case 'popupmenu'
        jComboBox = findjobj_fast(control);
        jComboBox.setRenderer(ListBoxRenderer);
end
