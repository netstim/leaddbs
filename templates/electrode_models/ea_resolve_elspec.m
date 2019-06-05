function varargout=ea_resolve_elspec(varargin)
% This simple function outputs a cellarray of available electrode specs if
% nargin==0 and exports the current electrode specification if varargin{1}
% is an options struct with options.elmodel defined as a valid name of electrode string.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

varargout{1}={'Medtronic 3389', 'Medtronic 3387', 'Medtronic 3391', ...
              'Boston Scientific Vercise', 'Boston Scientific Vercise Directed', ...
              'St. Jude ActiveTip (6146-6149)','St. Jude ActiveTip (6142-6145)',...
              'St. Jude Directed 6172 (short)','St. Jude Directed 6173 (long)', ...
              'PINS L301', 'PINS L302', 'PINS L303', ...
              'SDE-08 S8 Legacy', 'SDE-08 S10 Legacy', 'SDE-08 S12 Legacy', 'SDE-08 S16 Legacy', ...
              'SDE-08 S8', 'SDE-08 S10', 'SDE-08 S12', 'SDE-08 S16', ...
              '2069-EPC-05C-35', '2069-EPC-15C-35', 'NeuroPace DL-344-3.5'};
varargout{2}={'medtronic_3389', 'medtronic_3387', 'medtronic_3391', ...
              'boston_vercise', 'boston_vercise_directed', ...
              'stjude_activetip_2mm','stjude_activetip_3mm', ...
              'stjude_directed_05','stjude_directed_15', ...
              'pins_l301', 'pins_l302', 'pins_l303', ...
              'sde_08_s8_legacy', 'sde_08_s10_legacy', 'sde_08_s12_legacy', 'sde_08_s16_legacy',...
              'sde_08_s8', 'sde_08_s10', 'sde_08_s12', 'sde_08_s16', ...
              'epc_05c', 'epc_15c', 'neuropace_dl_344_35'};

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
            elspec.tip_length=1.1;
            elspec.contact_spacing=0.5;
            elspec.numel=8;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10','K11','K12','K13','K14','K15'};
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
        case 'St. Jude Directed 6172 (short)'
            elspec.matfname='stjude_directed_05';
            elspec.lead_diameter=1.27;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.27;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.27;
            elspec.tip_color=0.3;
            elspec.tip_length=1;
            elspec.contact_spacing=0.5;
            elspec.numel=8;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10','K11','K12','K13','K14','K15'};
        case 'St. Jude Directed 6173 (long)'
            elspec.matfname='stjude_directed_15';
            elspec.lead_diameter=1.27;
            elspec.lead_color=0.7;
            elspec.contact_length=1.5;
            elspec.contact_diameter=1.27;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.27;
            elspec.tip_color=0.3;
            elspec.tip_length=1;
            elspec.contact_spacing=1.5;
            elspec.numel=8;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10','K11','K12','K13','K14','K15'};
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
        case 'SDE-08 S8 Legacy'
            elspec.matfname='sde_08_s8_legacy';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=1.5;
            elspec.numel=8;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7',...
                                 'K8','K9','K10','K11','K12','K13','K14','K15'};
        case 'SDE-08 S10 Legacy'
            elspec.matfname='sde_08_s10_legacy';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=1.5;
            elspec.numel=10;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9',...
                                 'K10','K11','K12','K13','K14','K15','K16','K17','K18','K19'};
        case 'SDE-08 S12 Legacy'
            elspec.matfname='sde_08_s12_legacy';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=1.5;
            elspec.numel=12;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10','K11',...
                                 'K12','K13','K14','K15','K16','K17','K18','K19','K20','K21','K22','K23'};
        case 'SDE-08 S16 Legacy'
            elspec.matfname='sde_08_s16_legacy';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.7;
            elspec.tip_length=1.5;
            elspec.contact_spacing=1.5;
            elspec.numel=16;
            elspec.tipiscontact=0;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10','K11','K12','K13','K14','K15',...
                                 'K16','K17','K18','K19','K20','K21','K22','K23','K24','K25','K26','K27','K28','K29','K30','K31'};
        case 'SDE-08 S8'
            elspec.matfname='sde_08_s8';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.3;
            elspec.tip_length=2;
            elspec.contact_spacing=1.5;
            elspec.numel=8;
            elspec.tipiscontact=1;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7',...
                                 'K8','K9','K10','K11','K12','K13','K14','K15'};
        case 'SDE-08 S10'
            elspec.matfname='sde_08_s10';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.3;
            elspec.tip_length=2;
            elspec.contact_spacing=1.5;
            elspec.numel=10;
            elspec.tipiscontact=1;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9',...
                                 'K10','K11','K12','K13','K14','K15','K16','K17','K18','K19'};
        case 'SDE-08 S12'
            elspec.matfname='sde_08_s12';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.3;
            elspec.tip_length=2;
            elspec.contact_spacing=1.5;
            elspec.numel=12;
            elspec.tipiscontact=1;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10','K11',...
                                 'K12','K13','K14','K15','K16','K17','K18','K19','K20','K21','K22','K23'};
        case 'SDE-08 S16'
            elspec.matfname='sde_08_s16';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.3;
            elspec.tip_length=2;
            elspec.contact_spacing=1.5;
            elspec.numel=16;
            elspec.tipiscontact=1;
            elspec.contactnames={'K0','K1','K2','K3','K4',...
                                 'K16','K17','K18','K19','K20','K21','K22','K23','K24','K25','K26','K27','K28','K29','K30','K31'};
        case '2069-EPC-05C-35'
            elspec.matfname='epc_05c';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.3;
            elspec.tip_length=2;
            elspec.contact_spacing=1.5;
            elspec.numel=5;
            elspec.tipiscontact=1;
            elspec.contactnames={'K0','K1','K2','K3','K4',...
                                 'K8','K9','K10','K11','K12'};
        case '2069-EPC-15C-35'
            elspec.matfname='epc_15c';
            elspec.lead_diameter=0.8;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=0.8;
            elspec.contact_color=0.3;
            elspec.tip_diameter=0.8;
            elspec.tip_color=0.3;
            elspec.tip_length=2;
            elspec.contact_spacing=1.5;
            elspec.numel=15;
            elspec.tipiscontact=1;
            elspec.contactnames={'K0','K1','K2','K3','K4','K5','K6','K7','K8','K9','K10','K11','K12','K13','K14',...
                                 'K16','K17','K18','K19','K20','K21','K22','K23','K24','K25','K26','K27','K28','K29','K30'};
        case 'NeuroPace DL-344-3.5'
            elspec.matfname='neuropace_dl_344_35';
            elspec.lead_diameter=1.27;
            elspec.lead_color=0.7;
            elspec.contact_length=2;
            elspec.contact_diameter=1.27;
            elspec.contact_color=0.3;
            elspec.tip_diameter=1.27;
            elspec.tip_color=0.7;
            elspec.tip_length=1.1;
            elspec.contact_spacing=1.5;
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
