function [pamlist, obj] = ea_discfibers_getpams_lateral(obj)
% Return list of VATs

% this is a special version where we convert regular bilateral input
% to unilateral N*2

% split patients into left and right
numPatient = length(obj.allpatients);
pamlist = cell(numPatient*2,1);  


disp('Construct PAM list...')

% IMPORTANT: if multiple pathways were used, fiberActivation files have been already merged
% in ea_discfibers_merge_pathways!
for sub=1:numPatient % Original VAT E-field

    [~,subj_tag,~] = fileparts(obj.M.patient.list{sub});
    subSimPrefix = [subj_tag, '_sim-'];
    

    % IMPORTANT: we are going to change the patient list, the original
    % one needs to be restored to recalc! This can be done by cutting off
    % the second half


    pamlist{sub,1} = [obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-R.mat'];


    obj.allpatients{sub+numPatient,1} = [obj.M.patient.list{sub},'_from_left'];
    obj.M.patient.list{sub+numPatient,1} = [obj.M.patient.list{sub},'_from_left'];
    % load left here for now, will be flipped later
    pamlist{sub+numPatient,1} = [obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-L.mat'];


end


% also update scores
% I will completely overwrite them here
obj.M.clinical.labels{1,1} = 'Brady_Unilateral';
obj.M.clinical.labels{1,2} = 'Rigidity_Unilateral';
obj.M.clinical.labels{1,3} = 'Tremor_Unilateral';
obj.M.clinical.labels{1,4} = 'Axial_Unilateral';  % does not make sense for these
obj.M.clinical.labels{1,5} = 'Gait_Unilateral';     % does not make sense for these

obj.M.patient.group = ones(numPatient*2,1);
obj.M.patient.group(numPatient+1:end,1) = 2;

% first left body side, then right
obj.M.clinical.vars{1,1} = [0.273
1.0
0.5
0.4
0.6
0.364
0.286
0.308
0.688
0.625
0.5
0.125
0.263
-0.143
0.583
0.286
0.286
0.286
0.286
0.636
0.167
0.833
0.333
0.182
0.583
0.556
0.375
0.333
0.263
0.6
0.4
-0.5
0.5
0.0
0.111
0.0
nan
0.7
0.417
0.733
0.5
0.0
0.467
0.533
0.875
1.0
0.429
0.0
0.0
0.4
0.364
0.286
-0.333
0.5
0.571
-0.1
1.0
0.667
0.5
0.429
0.625
0.385
0.563
0.417
0.727
0.286
-0.286
0.5
0.1
0.0];

obj.M.clinical.vars{1,2} =[0.0
nan
0.842
0.667
0.818
0.5
1.0
0.6
0.8
1.0
nan
0.313
0.438
0.333
1.0
0.5
1.0
1.0
1.0
1.0
0.8
0.6
0.6
0.75
1.0
0.667
0.4
0.8
0.6
0.615
1.0
0.0
0.778
0.5
0.625
0.333
1.0
0.706
0.923
1.0
0.583
0.333
0.4
0.5
0.6
1.0
0.417
0.278
0.143
0.778
0.0
0.0
0.6
0.667
1.0
0.5
0.333
0.867
0.8
0.0
1.0
0.0
0.857
0.4
0.471
0.75
0.0
0.714
1.0
0.5];

obj.M.clinical.vars{1,3} = [0.167
1.0
0.571
1.0
nan
nan
1.0
nan
1.0
1.0
0.353
0.724
0.0
nan
0.0
1.0
0.35
0.7
1.0
1.0
1.0
1.0
1.0
1.0
0.5
1.0
1.0
1.0
0.556
nan
nan
nan
0.944
0.0
1.0
0.0
nan
0.667
0.733
nan
nan
1.0
nan
0.714
1.0
0.783
0.478
-0.667
nan
0.0
0.0
0.167
0.636
1.0
1.0
1.0
1.0
1.0
1.0
0.625
1.0
1.0
1.0
0.333
nan
nan
nan
0.85
0.0
1.0];

% we should also update selection based on unilateral baseline

%obj.responsevarlabel = obj.M.clinical.labels{1,1};
%obj.responsevar = obj.M.clinical.vars{1,1};

end