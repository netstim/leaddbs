function feval_onButtonDownExact(callback,hObject)

assert(~iscell(hObject), 'possible varargin{:} bug');
hFigure = getParentFigure(hObject);
hObject_pos = getPositionOnFigure(hObject,'pixels');


selectiontype = get(hFigure,'selectiontype');
if strcmp( selectiontype , 'alt')   || strcmp( selectiontype , 'open')
    pos=getCurrentPoint(hFigure,'pixels');
    
    if all(  pos>=hObject_pos(1:2) & pos<=hObject_pos(3:4)+hObject_pos(1:2)  )
        
        feval(callback,hObject);
        
    end
end
end