function ea_settransparency(resultfig,togglestates)

xsliceplot=getappdata(resultfig,'xsliceplot');
ysliceplot=getappdata(resultfig,'ysliceplot');
zsliceplot=getappdata(resultfig,'zsliceplot');

if ~isempty(xsliceplot)
    if isvalid(xsliceplot)
    set(xsliceplot,'FaceAlpha',togglestates.xyztransparencies(1))
    end
end
if ~isempty(ysliceplot)
    if isvalid(ysliceplot)
    set(ysliceplot,'FaceAlpha',togglestates.xyztransparencies(2))
    end
end
if ~isempty(zsliceplot)
    if isvalid(zsliceplot)
    set(zsliceplot,'FaceAlpha',togglestates.xyztransparencies(3))
    end
end

end
