function  showcontextmenu( hObject,uicontextmenu )
%SHOWCONTEXTMENU 
% 
% showcontextmenu( hObject ) shows uicontextmenu assigned to hObject
% 
% showcontextmenu( hObject,uicontextmenu ) shows custom uicontextmenu
% in a drop-down-menu like way
% 



if ~exist('uicontextmenu','var'); uicontextmenu=get(hObject,'uicontextmenu');end;

if isempty(uicontextmenu)
    error('showcontextmenu:uicontextmenu','no uicontextmenu found');
end

if  isprop(hObject,'Units') %buttons, etc.
    hObject_pos = getPositionOnFigure(hObject,'pixels');  % drop-down-like menu for button
else
    switch get(hObject,'type')
        case {'uipushtool','uitoggletool'} % thanks to Jesse <biodafes@gmail.com> for most of this code
            toolbar_handle = get(hObject, 'parent');
            sibling_handles = getAllChildren(toolbar_handle); %uipushtools of built-in toolbars got 'HandleVisibility' off
            sibling_handles = flipud(sibling_handles);% strangely the order is reverse that of the toolbar   
            sibling_gotsep_txt=get(sibling_handles,'Separator');  %separator takes some space too, result is {'off','on'..}
            sibling_gotsep = cellfun(@(s) strcmp('on',s) ,  sibling_gotsep_txt); %same in logicals
            
            i = find(sibling_handles == hObject, 1);
            x = (i - 1) * 23 ... space taken by icons, not sure if this doesn't depend on DPI , platform, etc...
                +sum(sibling_gotsep(1:i)) * 8 ... space taken by separators
                + 1;% Jesse got this equation by eyeballing the results and knowing that the button was somewhere around 20 pixels wide
            %note that this would be inaccurate if toolbar's length exeeds figure width
            
            
            % toolbar_handle parent must be figure
            fig = get(toolbar_handle, 'parent');
            xywh = get(fig, 'position');
            
            % y is the height of the figure
            hObject_pos = [x xywh(4)];
        otherwise  %line and etc.
            hObject_pos = get(getParentFigure(hObject),'currentpoint'); %on click position
    end
end
pos = hObject_pos(1:2);
set(uicontextmenu,'Position',pos);
set(uicontextmenu,'Visible','on');
end

