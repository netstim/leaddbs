function success = ea_runcleartune(ipStruct,tractset)


patlist = ipStruct.patlist;
app.tractset = tractset;
app.tractset.CleartuneOptim = 1;
tractsetclone=ea_disctract;
app.tractset.copyobj(tractsetclone);
app.inputVars.MinCylindricCurr = ipStruct.MinCylindricCurr;
app.inputVars.MaxCylindricCurr = ipStruct.MaxCylindricCurr;
app.inputVars.MinSegmentCurr = ipStruct.MinSegmentCurr;
app.inputVars.MaxSegmentCurr = ipStruct.MaxSegmentCurr;
app.inputVars.startContact = ipStruct.Startatcontact;
app.Applypenaltyusing = 'Quadratic curve';
app.SweetspotamplitudeEditField.Value = ipStruct.SweetspotAmplitude;
app.symptomWeightVar = ipStruct.scores;
options = ea_setopts_local;
options.native = 0;
options.groupmode = 1;
options.groupid = 'cleartune';
options = ea_getptopts(patlist{1}, options);
[coords_mm]=ea_load_reconstruction(options);
app.inputVars.numContacts = size(coords_mm{1},1);
app.inputVars.modelVTA = ipStruct.StimulationVolumeModel;
modeFlag = ipStruct.modeFlag; %choice of 'none' or 'restricted'
levelR = ipStruct.levelR; %this is specific for each patient.
levelL = ipStruct.levelL;
if strcmp(ipStruct.ConstCurr,'mA')
    app.inputVars.constcurr = 1; %1 for mA, 0 for V
else
    app.inputVars.constcurr = 0;
end
app.PenaltyVal = ipStruct.SweetspotAmplitude;
for i=1:length(patlist)
    [startptsR,startptsL,lbR,lbL,ubR,ubL] = setBounds(modeFlag,patlist{i},ipStruct.levelR,ipStruct.levelL);
    app.startptsR = startptsR;
    app.startptsL = startptsL;
    app.lbR = lbR;
    app.lbL = lbL;
    app.ubR = ubR;
    app.ubL = ubL;
end
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
function [startptsR,startptsL,lbR,lbL,ubR,ubL] = setBounds(modeFlag,sub,levelR,levelL)
    reconstruction_file = dir([sub,filesep,'reconstruction',filesep,'sub-*desc-reconstruction.mat']);
    reconstruction_file_path = [reconstruction_file.folder,filesep,reconstruction_file.name];
    [reconst, ~, ~, ~] = ea_get_reconstruction(reconstruction_file_path);
    % do not update S here, just get the bounds in mA!
    [min_bound_per_contactR, max_bound_per_contactR, ~] = ea_get_currents_per_contact(app.inputVars.MinCylindricCurr,app.inputVars.MaxCylindricCurr, app.inputVars.MinSegmentCurr, app.inputVars.MaxSegmentCurr, reconst, 0, 0);
    min_bound_per_contactL = min_bound_per_contactR;
    max_bound_per_contactL = max_bound_per_contactR;    
    if strcmpi(modeFlag,'restricted')
        firstLevel = [1,2,3,4,8];
        secondLevel = [1,5,6,7,8];
        if levelR == 1
            min_bound_per_contactR(secondLevel) = 0;
            max_bound_per_contactR(secondLevel) = 0;
            startcontactR = 3;
        elseif levelR == 2
            min_bound_per_contactR(firstLevel) = 0;
            max_bound_per_contactR(firstLevel) = 0;
            startcontactR = 5;
        end
        if levelL == 1
            min_bound_per_contactL(secondLevel) = 0;
            max_bound_per_contactL(secondLevel) = 0;
            startcontactL = 3;
        elseif levelL == 2
            min_bound_per_contactL(firstLevel) = 0;
            max_bound_per_contactL(firstLevel) = 0;
            startcontactL = 5;
        end
       
    end
    startptsR = zeros(1,app.inputVars.numContacts);
    startptsL = zeros(1,app.inputVars.numContacts);
    
    if abs(max_bound_per_contactR(startcontactR)) > abs(min_bound_per_contactR(startcontactR))
        startptsR(startcontactR) = max_bound_per_contactR(startcontactR) / 1.0;
    else
        startptsR(startcontactR) = min_bound_per_contactR(startcontactR) / 1.0;
    end
    if abs(max_bound_per_contactL(startcontactL)) > abs(min_bound_per_contactL(startcontactL))
        startptsL(startcontactL) = max_bound_per_contactL(startcontactL) / 1.0;
    else
        startptsL(startcontactL) = min_bound_per_contactL(startcontactL) / 1.0;
    end
    lbR = min_bound_per_contactR;
    lbL = min_bound_per_contactL;
    ubR = max_bound_per_contactR;
    ubL = max_bound_per_contactL;
end
end