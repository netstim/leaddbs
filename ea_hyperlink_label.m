function ea_hyperlink_label(label, url, position)
% Creates a clickable hyperlink label
%
% USAGE:
%
%    ea_hyperlink_label(label, url, position)
%
% INPUTS:
%    label:         label text
%    url:           linked URL
%    position:      position of the label

labelStr = ['<html><a href="">' label '</a></html>'];
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
[hjLabel,~] = ea_javacomponent(jLabel, position, gcf);
bgcolor = num2cell(get(gcf, 'Color'));
hjLabel.setBackground(java.awt.Color(bgcolor{:}));
hjLabel.setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR));
hjLabel.setToolTipText(['Click to visit ' url]);
set(hjLabel, 'MouseClickedCallback', @(h,e)web(url, '-browser'))
