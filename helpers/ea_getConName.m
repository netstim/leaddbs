function [conName, isDirectional] = ea_getConName(elmodel, side, opts)
% Get contact name and directionality
arguments
    elmodel % Lead model name/options/elspec
    side % 1/2/R/RH/L/LH
    opts.showSideStr  {mustBeNumericOrLogical} = false
end

if ischar(elmodel) % Lead model specified
    options.elmodel = elmodel;
    options = ea_resolve_elspec(options);
    conNum = options.elspec.numContacts;
elseif isstruct(elmodel)
    if isfield(elmodel, 'elspec') % 'options' specified
        options = elmodel;
    elseif isfield(elmodel, 'numContacts') %  'elspec' specified
        options.elmodel = elmodel.elmodel;
        options.elspec = elmodel;
    end
    elmodel = options.elspec.elmodel;
    conNum = options.elspec.numContacts;
end

switch side
    case {1, 'R', 'RH'} % RH
        conName = options.elspec.contactnames(1:conNum);
    case {2, 'L', 'LH'} % LH
        conName = options.elspec.contactnames(conNum+1:end);
    otherwise
        conName = options.elspec.contactnames(1:conNum);
        conName = erase(conName, " (R)" | " (L)");
end

if ~opts.showSideStr
    conName = erase(conName, " (R)" | " (L)");
end

switch elmodel
    case {'Medtronic B33005', 'Medtronic B33015', ...
          'Boston Scientific Vercise Directed', ...
          'Abbott Directed 6172 (short)', 'Abbott Directed 6173 (long)'}
        isDirectional = ones(1, conNum);
        isDirectional([1,end]) = 0;

    case 'Boston Scientific Vercise Cartesia HX'
        isDirectional = ones(1, conNum);
        isDirectional(end-3:end) = 0;

    case 'Boston Scientific Vercise Cartesia X'
        isDirectional = ones(1, conNum);
        isDirectional(end) = 0;

    case 'Aleva directSTIM Directed'
        isDirectional = ones(1, conNum);

    otherwise
        isDirectional = zeros(1, conNum);
end
