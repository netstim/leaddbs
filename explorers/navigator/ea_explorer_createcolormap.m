function fibcmap = ea_explorer_createcolormap(obj,gradientLevel)
cmapShiftRatio = 0;
shiftedCmapStart = round(gradientLevel*cmapShiftRatio)+1;
shiftedCmapEnd = gradientLevel-round(gradientLevel*cmapShiftRatio);
shiftedCmapLeftEnd = gradientLevel/2-round(gradientLevel/2*cmapShiftRatio);
shiftedCmapRightStart = round(gradientLevel/2*cmapShiftRatio)+1;

if obj.thresholding.posvisible && obj.thresholding.negvisible
    cmap = ea_colorgradient(gradientLevel/2, obj.thresholding.negcolor, [1,1,1]);
    cmapLeft = ea_colorgradient(gradientLevel/2, obj.thresholding.negcolor, cmap(shiftedCmapLeftEnd,:));
    cmap = ea_colorgradient(gradientLevel/2, [1,1,1], obj.thresholding.poscolor);
    cmapRight = ea_colorgradient(gradientLevel/2, cmap(shiftedCmapRightStart,:), obj.thresholding.poscolor);
    fibcmap = [cmapLeft;cmapRight];  
elseif obj.thresholding.posvisible
    cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.thresholding.poscolor);
    fibcmap = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.thresholding.poscolor);  
elseif obj.thresholding.negvisible
    cmap = ea_colorgradient(gradientLevel, obj.thresholding.negcolor, [1,1,1]);
    fibcmap = ea_colorgradient(gradientLevel, obj.thresholding.negcolor, cmap(shiftedCmapEnd,:)); 
end
end
