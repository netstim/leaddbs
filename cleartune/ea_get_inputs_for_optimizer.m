function [inputs,ampl,perc_val] = ea_get_inputs_for_optimizer(patient_folder, XOptim, modelVTA,writeVTA,side)

% New way to defide the stimulation window
% just because all cathodic

if any(XOptim(:) > 0.0)  && ~strcmp(modelVTA,'Simbio/FieldTrip (Horn 2017)') && ~strcmp(modelVTA,'FieldTrip')
    disp("Only cathodic stimulation is supported!")
    return
end

cathode_currents = sum(XOptim(XOptim < 0.0)); 
anode_currents = sum(XOptim(XOptim > 0.0)); 

% the difference will be assigned to casing
ampl = max([abs(cathode_currents), anode_currents]);

perc_val = zeros(1,8);
% keep original sign for perc here
perc_val(1:size(XOptim,2)) = 100.0 * XOptim./ampl;

constcurr = 0;  % 0 - CC, 1 - VC (Lead-DBS notation)
inputs = {patient_folder,ampl,perc_val,constcurr,side,writeVTA,modelVTA};
end