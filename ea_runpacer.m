function [coords_mm,trajectory,markers]=ea_runpacer(options)

directory = [options.root,options.patientname,filesep];
ctnii=options.prefs.ctnii_coregistered; tmat=eye(4);
if exist([options.root,options.patientname,filesep,stripext(options.prefs.rawctnii_unnormalized),'2',stripext(options.prefs.prenii_unnormalized),'_ants1.mat'],'file'); % use unresliced version and apply matrix in RAM
    load([options.root,options.patientname,filesep,stripext(options.prefs.rawctnii_unnormalized),'2',stripext(options.prefs.prenii_unnormalized),'_ants1.mat'])
    tmat=ea_antsmat2mat(AffineTransform_float_3_3,fixed);
    ctnii=options.prefs.rawctnii_unnormalized;
end
elecmodels=PaCER([options.root,options.patientname,filesep,ctnii],'finalDegree',1,'electrodeType',ea_mod2pacermod(options.elmodel));

for side=options.sides
    coords_mm{side}=[tmat*[elecmodels{side}.getContactPositions3D,ones(size(elecmodels{side}.getContactPositions3D,1),1)]']';
    coords_mm{side}=coords_mm{side}(:,1:3);
    for dim=1:3
        trajectory{side}(:,dim)=linspace(coords_mm{side}(1,dim),coords_mm{side}(1,dim)+10*(coords_mm{side}(1,dim)-coords_mm{side}(end,dim)),20);
    end
    
    markers(side).head=coords_mm{side}(1,:);
    markers(side).tail=coords_mm{side}(4,:);
    normtrajvector{side}=(coords_mm{side}(1,:)-coords_mm{side}(end,:))/...
        norm((coords_mm{side}(1,:)-coords_mm{side}(end,:)));
    orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);
    
    markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
    markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality
end


ea_methods(options,...
    ['DBS-Electrodes were automatically pre-localized in native & template space using the PaCER algorithm',...
    ' (Husch et al., 2017; http://adhusch.github.io/PaCER/).'],...
    {'Husch, A., Petersen, M. V., Gemmar, P., Goncalves, J., & Hertel, F. (2017). PaCER - A fully automated method for electrode trajectory and contact reconstruction in deep brain stimulation. NeuroImage. Clinical, 17, 80?89. http://doi.org/10.1016/j.nicl.2017.10.004'});

function fn=stripext(fn)
[~,fn]=fileparts(fn);


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
