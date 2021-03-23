function COG = ea_diode_calculateCOG(data,Xslice,Yslice,Zslice)  
    data(isnan(data)) = 0;
    values = data(:)';
%     [x,y] = ind2sub(size(data),[1:numel(data)]);
%     xval = Xslice(x,y).*values;
%     yval = Yslice(x,y).*values;
%     zval = Zslice(x,y).*values;

    
    xval = Xslice(:)'.*values;
    yval = Yslice(:)'.*values;
    zval = Zslice(:)'.*values;
    
    COG = [sum(xval) ./ sum(values); sum(yval) ./ sum(values);  sum(zval) ./ sum(values)];
end