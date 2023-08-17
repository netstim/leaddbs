function [coords_mm,trajectory,markers] = ea_runpacer(options)

[tmat,ctnii] = ea_getrawct2preniimat(options);
niiCTSPM = NiftiModSPM(ctnii); % load nifti using SPM instead of PaCER default Nifti Toolbox

if ~isfile(options.subj.recon.rawCTMask)
    ea_genctmask(options);
end
switch options.prefs.reco.mancoruse
    case 'postop'
        brainMask = options.subj.recon.rawCTMask;
    case 'rpostop'
        brainMask = options.subj.recon.anchorNativeMask;
end

elecmodels=PaCER(niiCTSPM,'finalDegree',1,'electrodeType',ea_mod2pacermod(options.elmodel), 'brainMask', brainMask);
disp('======== PaCER reconstruction finished. Converting PaCER reconstructions to LeadDBS. ========')

if(length(elecmodels) ~= length(options.sides))
   error(['PaCER returned a different number of electrodes than expected by LeadDBS! ' ...
       'In most cases this indicates an error in PaCER preprocessing (brain mask estimation) ' ...
       'due to untypical CT data. Please provide a brain mask to PaCER in this case using the mask parameter.']);
end

for iside=1:length(options.sides)
    side=options.sides(iside);

    coords_mm{side}=[tmat*[elecmodels{iside}.getContactPositions3D,ones(size(elecmodels{iside}.getContactPositions3D,1),1)]']';%#ok<NBRAK,AGROW>
    coords_mm{side}=coords_mm{side}(:,1:3);%#ok<AGROW>
    for dim=1:3
        trajectory{side}(:,dim)=linspace(coords_mm{side}(1,dim),coords_mm{side}(1,dim)+10*(coords_mm{side}(1,dim)-coords_mm{side}(end,dim)),20);%#ok<AGROW>
    end

    markers(side).head=coords_mm{side}(1,:); %#ok<AGROW>
    markers(side).tail=coords_mm{side}(4,:); %#ok<AGROW>
    [xunitv, yunitv] = ea_calcxy(coords_mm{side}(1,:), coords_mm{side}(end,:));
    markers(side).x=coords_mm{side}(1,:)+xunitv*(options.elspec.lead_diameter/2);
    markers(side).y=coords_mm{side}(1,:)+yunitv*(options.elspec.lead_diameter/2);
end


ea_methods(options,...
    ['DBS-Electrodes were automatically pre-localized in native & template space using the PaCER algorithm',...
    ' (Husch et al., 2017; http://adhusch.github.io/PaCER/).'],...
    {'Husch, A., Petersen, M. V., Gemmar, P., Goncalves, J., & Hertel, F. (2017). PaCER - A fully automated method for electrode trajectory and contact reconstruction in deep brain stimulation. NeuroImage. Clinical, 17, 80?89. http://doi.org/10.1016/j.nicl.2017.10.004'});


function model=ea_mod2pacermod(model)
% current dictionary to translate between Lead-DBS and PaCER nomenclature.
% Hoping to standardize this in the future.
switch model
    case 'Medtronic 3389'
        % pass through (same nomenclature)
    case 'Medtronic 3387'
        % pass through (same nomenclature)
    case 'Boston Scientific Vercise Directed'
        model='Boston Vercise Directional';
    otherwise
        model=''; % 'Unkown Electrode Type'
end
