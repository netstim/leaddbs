function ea_discfibers_showroi(obj)
% Function to draw seeds/roi/VTAs (from E-fields) in Lead-DBS
% fiberfiltering explorer



if ~obj.roivisible
    if ~isempty(obj.roidata)
        try
            for side=1:length(obj.roidata.nimage)
                try delete(obj.roidrawobject{side}.toggleH); end
                try delete(obj.roidrawobject{side}); end
            end
        end
    end
    obj.roiprotocol.drawn=0;
    return
end
if ea_ddf_sr_need_redraw(obj)
    prefs=ea_prefs;
    if ea_ddf_sr_need_recalcN(obj)
        ea_discfibers_roi_nimage_sel(obj); % run imcalc to generate N image from selection
    end

    % show ROI / N-map:
    pobj.plotFigureH = obj.resultfig;
    addht = getappdata(obj.resultfig,'addht');
    if isempty(addht)
        addht = uitoolbar(obj.resultfig);
        setappdata(obj.resultfig, 'addht', addht);
    end

    for side=1:length(obj.roidata.nimage)
        pobj.htH = addht;
        pobj.color = prefs.fibfilt.roi.color;
        pobj.alpha = prefs.fibfilt.roi.alpha;
        pobj.nii = obj.roidata.nimage_sel{side};
        pobj.threshold=0.5;
        try delete(obj.roidrawobject{side}); end
        obj.roidrawobject{side}=ea_roi('N-Image', pobj);
    end
end

% store protocol to check whether we need to redraw next time
obj.roiprotocol.patsel=obj.patientselection;
obj.roiprotocol.drawn=1;


function redraw=ea_ddf_sr_need_redraw(obj)
redraw=1;
if isfield(obj.roiprotocol, 'patsel')
    if isequal(obj.roiprotocol.patsel,obj.patientselection) && obj.roiprotocol.drawn==1
        redraw=0;
    end
end


function recalcN=ea_ddf_sr_need_recalcN(obj)
recalcN=1;
if isfield(obj.roiprotocol, 'patsel')
    if isequal(obj.roiprotocol.patsel,obj.patientselection)
        recalcN=0;
    end
end
