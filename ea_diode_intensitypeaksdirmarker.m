function [sumintensity] = ea_diode_intensitypeaksdirmarker(intensity,angles)
    %% this function detects 'noPeaks' number of intensity peaks. peaks are constrained to be at 360°/noPeaks angles to each other.
    %% Function runs a noPeaks * (360°/noPeaks) array over the intensity-profile and finds the angle at which the sum of all peaks is highest.
    peak = round(rad2deg(angles) +1);
    peak(find(peak<1)) = peak(find(peak<1)) +360;
    peak(find(peak>360)) = peak(find(peak>360)) -360;    
    sumintensity = sum(intensity(peak));    
end