function new_tractset = ea_unifiedmapping_update_settings(explorer,new_settings)
    % Updates the given tractset object based on the
    % struct passed into the function

    explorer.calcsettings.calcthreshold = new_settings.calcthreshold;
    explorer.posvisible = new_settings.posvisible;
    explorer.negvisible = new_settings.negvisible;
    explorer.showposamount = new_settings.showposamount;
    explorer.shownegamount = new_settings.shownegamount;
    explorer.statsettings.connthreshold = new_settings.connthreshold;
    explorer.statsettings.efieldthreshold = new_settings.efieldthreshold;
    explorer.statsettings.statmetric = new_settings.statmetric;
    if isfield(new_settings, 'stimulationmodel')
        explorer.statsettings.stimulationmodel = new_settings.stimulationmodel;
    end
    explorer.corrtype = new_settings.corrtype;
    explorer.statsettings.efieldmetric = new_settings.efieldmetric;
    explorer.multitractmode = new_settings.multitractmode;
    explorer.numpcs = new_settings.numpcs;
    explorer.doactualprediction = new_settings.doactualprediction;
    explorer.predictionmodel = new_settings.predictionmodel;
    explorer.showsignificantonly = new_settings.showsignificantonly;
    explorer.alphalevel = new_settings.alphalevel;
    explorer.multcompstrategy = new_settings.multcompstrategy;
    explorer.basepredictionon = new_settings.basepredictionon;
    explorer.mirrorsides = new_settings.mirrorsides;
    explorer.modelNormalization = new_settings.modelNormalization;
    if isfield(new_settings, 'mask_vta_fp')
        explorer.mask_vta_fp = new_settings.mask_vta_fp;
    end
    explorer.numBins = new_settings.numBins;
    explorer.Nperm = new_settings.Nperm;
    explorer.kfold = new_settings.kfold;
    explorer.Nsets = new_settings.Nsets;
    explorer.adjustforgroups = new_settings.adjustforgroups;
    if isfield(new_settings, 'connectivity_type')
        explorer.calcsettings.connectivity_type = new_settings.connectivity_type;
    end
    if isfield(new_settings, 'switch_connectivity')
        explorer.calcsettings.switch_connectivity = new_settings.switch_connectivity;
    end
    if isfield(new_settings, 'calcmethod')
        explorer.calcsettings.calcmethod = new_settings.calcmethod;
    end
    if isfield(new_settings, 'calcspace')
        explorer.calcsettings.calcspace = new_settings.calcspace;
    end

    new_tractset = explorer;
end
