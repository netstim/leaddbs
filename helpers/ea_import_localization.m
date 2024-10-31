function ea_import_localization(dataset, subjID, coords, space, elmodel)
% Import localization from other tool into Lead-DBS subject folder.
arguments
    dataset {mustBeTextScalar} % Dataset path
    subjID  {mustBeTextScalar} % Subject ID
    coords
    space {mustBeMember(space, {'native', 'mni'})} = 'mni' % Space of the coordinates
    elmodel {mustBeTextScalar} = '' % Electrode model
end

models = ea_resolve_elspec;
if isempty(elmodel)
    index = listdlg('PromptString', 'Which electrode model it is?', 'ListString', models, 'SelectionMode', 'single', 'CancelString', 'Cancel');
     if ~isempty(index)
        elmodel = models{index};
    else
        return;
    end
elseif ~ismember(elmodel, models)
    ea_error('Specified model is not supported in Lead-DBS!', simpleStack=true);
end

options.elmodel = elmodel;
options = ea_resolve_elspec(options);

if isnumeric(coords)
    coords = {coords};
end

for i=1:length(coords)
    reco.props(i).elmodel = elmodel;
    reco.props(i).manually_corrected = 1;
    reco.(space).coords_mm{i} = coords{i};
    reco.(space).markers(i).head = coords{i}(1, :);
    switch elmodel
        case {'Medtronic B33005'
              'Medtronic B33015'
              'Boston Scientific Vercise Directed'
              'Abbott Directed 6172 (short)'
              'Abbott Directed 6173 (long)'}
            reco.(space).markers(i).tail = coords{i}(8,:);
        case {'Boston Scientific Vercise Cartesia HX'
              'Boston Scientific Vercise Cartesia X'}
            reco.(space).markers(i).head = mean(coords{i}(1:3,:));
            reco.(space).markers(i).tail = mean(coords{i}(10:12,:));
        otherwise
            reco.(space).markers(i).tail = coords{i}(4,:);
    end


    [xunitv, yunitv] = ea_calcxy(reco.(space).markers(i).head, reco.(space).markers(i).tail);
    reco.(space).markers(i).x = reco.(space).markers(i).head + xunitv*(options.elspec.lead_diameter/2);
    reco.(space).markers(i).y = reco.(space).markers(i).head + yunitv*(options.elspec.lead_diameter/2);
end

[~, reco.(space).trajectory] = ea_resolvecoords(reco.(space).markers, elmodel);

bids = BIDSFetcher(dataset);
subj = bids.getSubj(erase(subjID, textBoundary("start") + 'sub-'));

ea_mkdir(fileparts(subj.recon.recon));
save(subj.recon.recon, 'reco');
