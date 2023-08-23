function success = ea_runcleartune(ipStruct)


patlist = ipStruct.patlist;
app.tractset = ipStruct.tractset;
app.tractset.CleartuneOptim = 1;
tractsetclone=ea_disctract;
app.tractset.copyobj(tractsetclone);
app.inputVars.MinCylindricCurr = ipStruct.MinCylindricCurr;
app.inputVars.MaxCylindricCurr = ipStruct.MaxCylindricCurr;
app.inputVars.MinSegmentCurr = ipStruct.MinSegmentCurr;
app.inputVars.MaxSegmentCurr = ipStruct.MaxSegmentCurr;
app.inputVars.startContact = ipStruct.Startatcontact;
app.Applypenaltyusing = 'Quadratic Curve';
app.SweetspotamplitudeEditField.Value = 2;
app.symptomWeightVar = ipStruct.scores;
options = ea_setopts_local;
options.native = 0;
options.groupmode = 1;
options.groupid = 'cleartune';
options = ea_getptopts(patlist{1}, options);
[coords_mm]=ea_load_reconstruction(options);
app.inputVars.numContacts = size(coords_mm{1},1);
app.inputVars.modelVTA = ipStruct.StimulationVolumeModel;
if strcmp(ipStruct.ConstCurr,'mA')
    app.inputVars.constcurr = 1; %1 for mA, 0 for V
else
    app.inputVars.constcurr = 0;
end
app.PenaltyVal = ipStruct.SweetspotAmplitude;
ea_cleartune_cluster(tractsetclone,patlist,app);

function options=ea_setopts_local
    options.earoot=ea_getearoot;
    options.verbose=3;
    options.sides=1:2; 
    options.fiberthresh=1;
    options.writeoutstats=1;
    options.writeoutpm = 0;
    options.colormap=jet;
    options.d3.write=1;
    options.d3.prolong_electrode=2;
    options.d3.writeatlases=1;
    options.macaquemodus=0;
end
end