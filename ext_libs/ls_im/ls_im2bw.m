function BW=ls_im2bw(img,t)
    switch class(img)
        case {'double', 'single'}
            BW = (img >= t);
        case {'uint8'}
            BW = (img >= 255*t);
        case {'uint16'}
            BW = (img >= 65535*t);
        otherwise
            error('unsupported image class');
    end
end