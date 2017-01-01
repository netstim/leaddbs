function varargout=ea_resolve_elspec(varargin)
% This simple function outputs a cellarray of available electrode specs if
% nargin==0 and exports the current electrode specification if varargin{1}
% is an options struct with options.elmodel defined as a valid name of electrode string.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


varargout{1}={'Medtronic 3389', 'Medtronic 3387', 'Medtronic 3391', ...
              'Boston Scientific Vercise', 'Boston Scientific Vercise Directed', ...
              'St. Jude ActiveTip (6146-6149)','St. Jude ActiveTip (6142-6145)' ...
              'PINS L301', 'PINS L302' 'PINS L303'};
varargout{2}={'medtronic_3389', 'medtronic_3387', 'medtronic_3391', ...
              'boston_vercise', 'boston_vercise_directed', ...
              'stjude_activetip_2mm','stjude_activetip_3mm' ...
              'pins_l301', 'pins_l302', 'pins_l303'};
if ~nargin
    return
else
   options=varargin{1};
end

try
    switch options.elmodel
        case 'Medtronic 3389'
            elspec.matfname='medtronic_3389';
            elspec.lead_diameter=1.27;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.27;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.27;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=0.5;
            elspec.numel=4;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
        case 'Medtronic 3387'
            elspec.matfname='medtronic_3387';
            elspec.lead_diameter=1.27;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.27;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.27;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=1.5;
            elspec.numel=4;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
        case 'Medtronic 3391'
            elspec.matfname='medtronic_3391';
            elspec.lead_diameter=1.27;
            elspec.lead_color=0.7;
            elspec.contact_length=3;
            elspec.contact_diameter=1.27;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.27;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=4;
            elspec.numel=4;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
        case 'Boston Scientific Vercise'
            elspec.matfname='boston_vercise';
            elspec.lead_diameter=1.3;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.3;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.3;
            elspec.tip_color=0.7;
            elspec.tip_length=1.3;
            elspec.contact_spacing=0.5;
            elspec.numel=8;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10','K11','K12','K13','K14','K15'};
        case 'St. Jude ActiveTip (6146-6149)'
            elspec.matfname='stjude_activetip_2mm';
            elspec.lead_diameter=1.27;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.27;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.27;
            elspec.tip_color=0.3;
            elspec.tip_length=1.5;
            elspec.contact_spacing=0.5;
            elspec.numel=4;
            elspec.tipiscontact=1;
            elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
        case 'St. Jude ActiveTip (6142-6145)'
            elspec.matfname='stjude_activetip_3mm';
            elspec.lead_diameter=1.27;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.27;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.27;
            elspec.tip_color=0.3;
            elspec.tip_length=1.5;
            elspec.contact_spacing=1.5;
            elspec.numel=4;
            elspec.tipiscontact=1;
            elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
        case 'Boston Scientific Vercise Directed'
            elspec.matfname='boston_vercise_directed';
            elspec.lead_diameter=1.3;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.3;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.3;
            elspec.tip_color=0.3;
            elspec.tip_length=1.5;
            elspec.contact_spacing=0.5;
            elspec.numel=8; % correct here since the directional leads will be inflated lateron.
            elspec.tipiscontact=1;
            elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
        case 'PINS L301'
            elspec.matfname='pins_l301';
            elspec.lead_diameter=1.3;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.3;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.3;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=0.5;
            elspec.numel=4;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
        case 'PINS L302'
            elspec.matfname='pins_l302';
            elspec.lead_diameter=1.3;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.3;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.3;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=1.5;
            elspec.numel=4;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
        case 'PINS L303'
            elspec.matfname='pins_l303';
            elspec.lead_diameter=1.3;
            elspec.lead_color=0.7;
            elspec.contact_length=3;
            elspec.contact_diameter=1.3;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.3;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=3;
            elspec.numel=4;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
    end
catch
%    warning('No electrode model specified. Using Medtronic 3389.');
    elspec.matfname='medtronic_3389';
    elspec.lead_diameter=1.27;
    elspec.lead_color=0.7;
    elspec.contact_length=1.5;
    elspec.contact_diameter=1.27;
    elspec.contact_color=0.3;
    elspec.tip_diameter=1.27;
    elspec.tip_color=0.7;
    elspec.tip_length=1.5;
    elspec.contact_spacing=0.5;
    elspec.numel=4;
    elspec.tipiscontact=0;
    elspec.contactnames={'K0','K1','K2','K3','K8','K9','K10','K11'};
end

if ~isfield(elspec,'eldist')
    elspec.eldist=elspec.contact_spacing+elspec.contact_length;
end

try
    options.elspec=elspec;
catch
    keyboard
end

varargout{1}=options;
