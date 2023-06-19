function [contactnames,directional] = ea_getelcontactnames(elspec,side)
% Get contact name and directionality 

switch side
    case 1 % RH
        contactnames = elspec.contactnames(1:length(elspec.contactnames)/2);
    case 2 % LH
        contactnames = elspec.contactnames((length(elspec.contactnames)/2)+1:end);
end

switch elspec.elmodel
    case {'Medtronic B33005', 'Medtronic B33015', ...
          'Boston Scientific Vercise Directed', ...
          'Abbott Directed 6172 (short)', 'Abbott Directed 6173 (long)'}
        directional = ones(1, elspec.numel);
        directional([1,end]) = 0;

    case 'Boston Scientific Vercise Cartesia HX'
        directional = ones(1, elspec.numel);
        directional(end-3:end) = 0;

    case 'Boston Scientific Vercise Cartesia X'
        directional = ones(1, elspec.numel);
        directional(end) = 0;

    case 'Aleva directSTIM Directed'
        directional = ones(1, elspec.numel);

    otherwise
        directional = zeros(1, elspec.numel);
end
