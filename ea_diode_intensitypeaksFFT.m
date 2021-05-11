function [peak,sprofil] = ea_diode_intensitypeaksFFT(intensity,noPeaks)
    %% this function detects 'noPeaks' number of intensity peaks. peaks are constrained to be at 360°/noPeaks angles to each other.
    %% Function runs a noPeaks * (360°/noPeaks) array over the intensity-profile and finds the angle at which the sum of all peaks is highest.
    fftint = fft(intensity);
    fftpart = fftint(noPeaks+1);
    amplitude = abs(fftpart);
    phase = -asin(real(fftpart) / amplitude);
    if imag(fftpart) > 0
        if real(fftpart) > 0
            phase = -pi -phase;
        else
            phase = pi -phase;
        end
    end
    clear amplitude
    amplitude = (max(intensity) + abs(min(intensity))) / 2;
    level = max(intensity) - amplitude;
    for k = 1:360
%        sprofil(k) =amplitude * sin(deg2rad(2*k)-phase) + level;
       sprofil(k) =amplitude * sin(deg2rad(noPeaks*k)-phase) + level;
    end
    
    for k = 1:noPeaks
        peak(k) = (k-1)*(360/noPeaks) +1;
    end
    
    for k = 1:(360/noPeaks)
        sumintensity(k) = sum(sprofil(peak));
        peak = peak +1;
    end
    
    [~,maxpeak] = max(sumintensity);
    
    for k = 1:noPeaks
        peak(k) = maxpeak + (k-1)*(360/noPeaks);
    end
end