function mycolor = ea_explorer_value2color(value,colormap,coloraxis)
        minval = coloraxis(1);
        maxval = coloraxis(2);
        diffval = coloraxis(2)-coloraxis(1);
        coloridx = round((value-minval)  .* ((size(colormap,1)-1)/diffval)+1);
        coloridx(coloridx < 1) = 1;
        coloridx(coloridx > size(colormap,1)) = size(colormap,1);
        coloridx(isnan(coloridx)) = 1;
        mycolor = colormap(coloridx,:);
        arraysize = size(value);
        mycolor = reshape(mycolor,[arraysize,3]);
end