function [Efield] = ea_get_efield(perc,standard_efield,amp,conductivity,amp_mode,impedence)
% inputs:
% percentage on each contact
% standard_efield: the standard efield is provided in the resources folder
% of the fastfield
% amp: stimulation amplitude
% Electrode_type: the name of the stimulation electrode
% conductivity: the brain conductivity

%ouput:
% Efield: the electric field with the defiend parameters

if strcmp(amp_mode,'V')
    new_amp = (amp/impedence)*1000;
    amp = new_amp;
    amp_mode = 'mA';
end

perc=perc/100;

% calculate the efield for this configuration
Efield = zeros(size(standard_efield{1}));

for len = 1:length(standard_efield)
    Efield = Efield + (perc(len)*standard_efield{len,1});
end

% scale the conductivity by 0.14 because the standard efield is generated
% with conductivity of 0.1
if strcmp(amp_mode,'mA')
    scale_con= 0.1 / conductivity;
    Efield = Efield*scale_con;
end

% scale the electric field by amplitude, the standard efield is genereated
% for amplitude of 1mA
Efield = Efield*amp;

Efield(isnan(Efield))=0;
Efield(Efield>10000) = 10000;

% Smooth data
Efield = smooth3(Efield,'box', 9);
