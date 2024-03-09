function adaptiveThreshold = ea_get_adaptiveEthreshold(pulse_width)
    
    % get a binarization threshold for E-field based on pulse width
    % the thresholded is described by a "strength-duration" curve
    % a / (bx + c) + d fit based on Astrom 2014 (10.1109/TBME.2014.2363494)
    % for 3 um fiber diameter 
    %  PWs:      30,   60,   90,  120
    %  ||E||: 376.0 240.0 185.0 157.0

    a = -59806189;
    b = -4810;
    c = -45282;
    d = 60.56;

    adaptiveThreshold = a / (b * pulse_width + c) + d;

end
