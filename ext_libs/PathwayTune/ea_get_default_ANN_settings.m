function ANN_default_settings = ea_get_default_ANN_settings()

    % pathway settings
    ANN_default_settings.connectome = 'Multi-Tract: M1LowerExDS';
    ANN_default_settings.pathway = 'M1_cf_lowerex';
    ANN_default_settings.fiberDiameter = 4.0;           % in um 
    ANN_default_settings.axonLength = 15.0;             % in mm
    ANN_default_settings.axonModel = 'McNeal1976';      % OSS-DBS axon models
    ANN_default_settings.optim_space = 'native';        % native or MNI

    % Current settings
    ANN_default_settings.MinCylindricCurrent = -7.0;    % in A (- for cathodes)
    ANN_default_settings.MaxCylindricCurrent = 0.0;
    ANN_default_settings.MinSegmentedCurrent = -7.0;
    ANN_default_settings.MaxSegmentedCurrent = 0.0;

    % ANN errors
    ANN_default_settings.ErrorThreshold = 10.0;                              % ANN error of percent activation (also used in the capsule model)
    ANN_default_settings.SideEffectErrorThreshold = 5.0;                     % ANN error of percent activation for side-effect pathways
    
    % Optimizer settings (not relevant for the capsule model)
    ANN_default_settings.NumberofIterations = 5000;                          % number of optimizer iterations
    ANN_default_settings.OptAlgorithm = 'Dual Annealing';  
    ANN_default_settings.DistanceMetric = 'Canberra'; 
end
