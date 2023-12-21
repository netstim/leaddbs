function new_tractset = ea_update_settings(tractset,new_settings)
    % Updates the given tractset object based on the
    % struct passed into the function

    tractset.calcthreshold = new_settings.calcthreshold;
    tractset.posvisible = new_settings.posvisible;
    tractset.negvisible = new_settings.negvisible;
    tractset.showposamount = new_settings.showposamount;
    tractset.shownegamount = new_settings.shownegamount;
    tractset.connthreshold = new_settings.connthreshold;
    tractset.efieldthreshold = new_settings.efieldthreshold;
    tractset.statmetric = new_settings.statmetric;
    tractset.corrtype = new_settings.corrtype;
    tractset.efieldmetric = new_settings.efieldmetric;
    tractset.multitractmode = new_settings.multitractmode;
    tractset.numpcs = new_settings.numpcs;
    tractset.doactualprediction = new_settings.doactualprediction;
    tractset.predictionmodel = new_settings.predictionmodel;
    tractset.showsignificantonly = new_settings.showsignificantonly;
    tractset.alphalevel = new_settings.alphalevel;
    tractset.multcompstrategy = new_settings.multcompstrategy;
    tractset.basepredictionon = new_settings.basepredictionon;
    tractset.mirrorsides = new_settings.mirrorsides;
    tractset.modelNormalization = new_settings.modelNormalization;
    tractset.numBins = new_settings.numBins;
    tractset.Nperm = new_settings.Nperm;
    tractset.kfold = new_settings.kfold;
    tractset.Nsets = new_settings.Nsets;
    tractset.adjustforgroups = new_settings.adjustforgroups;

    new_tractset = tractset;
end