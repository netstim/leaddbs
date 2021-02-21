function varargout=ea_resolve_elspec(varargin)
% This simple function outputs a cellarray of available electrode specs if
% nargin==0 and exports the current electrode specification if varargin{1}
% is an options struct with options.elmodel defined as a valid name of electrode string.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ~nargin
    varargout{1}={'Medtronic 3389', 'Medtronic 3387', 'Medtronic 3391', ...
        'Boston Scientific Vercise', 'Boston Scientific Vercise Directed', ...
        'St. Jude ActiveTip (6146-6149)','St. Jude ActiveTip (6142-6145)', ...
        'St. Jude Directed 6172 (short)','St. Jude Directed 6173 (long)', ...
        'PINS Medical L301', 'PINS Medical L302', 'PINS Medical L303', .....
        'SDE-08 S8 Legacy', 'SDE-08 S10 Legacy', 'SDE-08 S12 Legacy', 'SDE-08 S16 Legacy', ...
        'SDE-08 S8', 'SDE-08 S10', 'SDE-08 S12', 'SDE-08 S14', 'SDE-08 S16', ...
        'PMT 2102-16-092', 'PMT 2102-16-093', 'PMT 2102-16-142', ...
        '2069-EPC-05C-35', '2069-EPC-15C-35', 'NeuroPace DL-344-3.5', 'NeuroPace DL-344-10', ...
        'DIXI D08-05AM', 'DIXI D08-08AM', 'DIXI D08-10AM', 'DIXI D08-12AM', 'DIXI D08-15AM', 'DIXI D08-18AM', ...
        'AdTech SD10R-SP05X Choi', 'AdTech RD10R-SP03X', 'AdTech BF08R-SP05X', 'AdTech BF08R-SP21X', 'AdTech BF08R-SP61X', ...
        'ELAINE Rat Electrode', 'NuMed Mini Lead'}';
    varargout{2}={'medtronic_3389', 'medtronic_3387', 'medtronic_3391', ...
        'boston_vercise', 'boston_vercise_directed', ...
        'stjude_activetip_2mm','stjude_activetip_3mm', ...
        'stjude_directed_05','stjude_directed_15', ...
        'pins_l301', 'pins_l302', 'pins_l303', ...
        'sde_08_s8_legacy', 'sde_08_s10_legacy', 'sde_08_s12_legacy', 'sde_08_s16_legacy',...
        'sde_08_s8', 'sde_08_s10', 'sde_08_s12', 'sde_08_s14', 'sde_08_s16', ...
        'pmt_2102_16_092', 'pmt_2102_16_093', 'pmt_2102_16_142', ...
        'epc_05c', 'epc_15c', 'neuropace_dl_344_35', ...
        'dixi_d08_05am', 'dixi_d08_08am', 'dixi_d08_10am', 'dixi_d08_12am', 'dixi_d08_15am', 'dixi_d08_18am', ...
        'adtech_sd10r_sp05x_choi', 'adtech_rd10r_sp03x', 'adtech_bf08r_sp05x', 'adtech_bf08r_sp21x', 'adtech_bf08r_sp61x', ...
        'elaine_rat_electrode', 'numed_minilead'}';
    return
else
    options=varargin{1};
end

if ~isfield(options, 'elmodel')
    try
        load([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'reco');
        elmodel = ea_get_first_notempty_elmodel(reco.props);
    catch
        %no model was found
        warning('No electrode model specified. Using Medtronic 3389.');
        elmodel = 'Medtronic 3389';
    end
else
    elmodel = options.elmodel;
end

switch elmodel
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
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
        elspec.contactnames={'K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)','K15 (R)','K16 (R)','K1 (L)','K2 (L)','K3 (L)','K4 (L)','K5 (L)','K6 (L)','K7 (L)','K8 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
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
        elspec.numel=8;
        elspec.tipiscontact=1;
        elspec.markerpos = 11;
        elspec.markerlen = 3;
        elspec.contactnames={'K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)','K15 (R)','K16 (R)','K1 (L)','K2 (L)','K3 (L)','K4 (L)','K5 (L)','K6 (L)','K7 (L)','K8 (L)'};
        elspec.isdirected=1;
        elspec.etagenames{1}={'K9 (R)','K10-12 (R)','K13-15 (R)','K16 (R)'};
        elspec.etagenames{2}={'K1 (L)','K2-4 (L)','K5-7 (L)','K8 (L)'};
        elspec.etageidx={1,2:4,5:7,8};
        elspec.forstimulation=1;
    case 'St. Jude ActiveTip (6146-6149)'
        elspec.matfname='stjude_activetip_2mm';
        elspec.lead_diameter=1.4;
        elspec.lead_color=0.7;
        elspec.contact_length=1.5;
        elspec.contact_diameter=1.4;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.4;
        elspec.tip_color=0.3;
        elspec.tip_length=3.0;
        elspec.contact_spacing=0.5;
        elspec.numel=4;
        elspec.tipiscontact=1;
        elspec.contactnames={'K1 (R)','K2 (R)','K3 (R)','K4 (R)','K1 (L)','K2 (L)','K3 (L)','K4 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'St. Jude ActiveTip (6142-6145)'
        elspec.matfname='stjude_activetip_3mm';
        elspec.lead_diameter=1.4;
        elspec.lead_color=0.7;
        elspec.contact_length=1.5;
        elspec.contact_diameter=1.4;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.4;
        elspec.tip_color=0.3;
        elspec.tip_length=3;
        elspec.contact_spacing=1.5;
        elspec.numel=4;
        elspec.tipiscontact=1;
        elspec.contactnames={'K1 (R)','K2 (R)','K3 (R)','K4 (R)','K1 (L)','K2 (L)','K3 (L)','K4 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
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
        elspec.markerpos = 10.75;
        elspec.markerlen = 1.5;
        elspec.contactnames={'K1 (R)','K2A (R)','K2B (R)','K2C (R)','K3A (R)','K3B (R)','K3C (R)','K4 (R)','K1 (L)','K2A (L)','K2B (L)','K2C (L)','K3A (L)','K3B (L)','K3C (L)','K4 (L)'};
        elspec.isdirected=1;
        elspec.etagenames{1}={'K1 (R)','K2 (R)','K3 (R)','K4 (R)'};
        elspec.etagenames{2}={'K1 (L)','K2 (L)','K3 (L)','K4 (L)'};
        elspec.etageidx={1,2:4,5:7,8};
        elspec.forstimulation=1;
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
        elspec.markerpos = 13.75;
        elspec.markerlen = 1.5;
        elspec.contactnames={'K1 (R)','K2A (R)','K2B (R)','K2C (R)','K3A (R)','K3B (R)','K3C (R)','K4 (R)','K1 (L)','K2A (L)','K2B (L)','K2C (L)','K3A (L)','K3B (L)','K3C (L)','K4 (L)'};
        elspec.isdirected=1;
        elspec.etagenames{1}={'K1 (R)','K2 (R)','K3 (R)','K4 (R)'};
        elspec.etagenames{2}={'K1 (L)','K2 (L)','K3 (L)','K4 (L)'};
        elspec.etageidx={1,2:4,5:7,8};
        elspec.forstimulation=1;
    case 'PINS Medical L301'
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'PINS Medical L302'
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'PINS Medical L303'
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)',...
            'K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)','K15 (R)',...
            'K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (L)','K27 (L)','K28 (L)','K29 (L)','K30 (L)','K31 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)',...
            'K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
    case 'SDE-08 S14'
        elspec.matfname='sde_08_s14';
        elspec.lead_diameter=0.8;
        elspec.lead_color=0.7;
        elspec.contact_length=2;
        elspec.contact_diameter=0.8;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.8;
        elspec.tip_color=0.3;
        elspec.tip_length=2;
        elspec.contact_spacing=1.5;
        elspec.numel=14;
        elspec.tipiscontact=1;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)',...
            'K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (L)','K27 (L)','K28 (L)','K29 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)','K15 (R)',...
            'K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (L)','K27 (L)','K28 (L)','K29 (L)','K30 (L)','K31 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
    case 'PMT 2102-16-092'
        elspec.matfname='pmt_2102_16_092';
        elspec.lead_diameter=0.8;
        elspec.lead_color=0.7;
        elspec.contact_length=2;
        elspec.contact_diameter=0.8;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.8;
        elspec.tip_color=0.3;
        elspec.tip_length=2;
        elspec.contact_spacing=1.97;
        elspec.numel=16;
        elspec.tipiscontact=1;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)','K15 (R)',...
            'K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (L)','K27 (L)','K28 (L)','K29 (L)','K30 (L)','K31 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
    case 'PMT 2102-16-093'
        elspec.matfname='pmt_2102_16_093';
        elspec.lead_diameter=0.8;
        elspec.lead_color=0.7;
        elspec.contact_length=2;
        elspec.contact_diameter=0.8;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.8;
        elspec.tip_color=0.3;
        elspec.tip_length=2;
        elspec.contact_spacing=2.43;
        elspec.numel=16;
        elspec.tipiscontact=1;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)','K15 (R)',...
            'K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (L)','K27 (L)','K28 (L)','K29 (L)','K30 (L)','K31 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
    case 'PMT 2102-16-142'
        elspec.matfname='pmt_2102_16_142';
        elspec.lead_diameter=0.8;
        elspec.lead_color=0.7;
        elspec.contact_length=2;
        elspec.contact_diameter=0.8;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.8;
        elspec.tip_color=0.3;
        elspec.tip_length=2;
        elspec.contact_spacing=3.53;
        elspec.numel=16;
        elspec.tipiscontact=1;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)','K15 (R)',...
            'K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (L)','K27 (L)','K28 (L)','K29 (L)','K30 (L)','K31 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)',...
            'K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (L)','K27 (L)','K28 (L)','K29 (L)','K30 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'NeuroPace DL-344-10'
        elspec.matfname='neuropace_dl_344_10';
        elspec.lead_diameter=1.27;
        elspec.lead_color=0.7;
        elspec.contact_length=2;
        elspec.contact_diameter=1.27;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.27;
        elspec.tip_color=0.7;
        elspec.tip_length=1.1;
        elspec.contact_spacing=8;
        elspec.numel=4;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'DIXI D08-05AM'
        elspec.matfname='dixi_d08_05am';
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)',...
            'K5 (L)','K6 (L)','K7 (L)','K8 (L)','K9 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'DIXI D08-08AM'
        elspec.matfname='dixi_d08_08am';
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'DIXI D08-10AM'
        elspec.matfname='dixi_d08_10am';
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
        elspec.isdirected=1;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
    case 'DIXI D08-12AM'
        elspec.matfname='dixi_d08_12am';
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)',...
            'K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'DIXI D08-15AM'
        elspec.matfname='dixi_d08_15am';
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
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)',...
            'K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (L)','K27 (L)','K28 (L)','K29 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'DIXI D08-18AM'
        elspec.matfname='dixi_d08_18am';
        elspec.lead_diameter=0.8;
        elspec.lead_color=0.7;
        elspec.contact_length=2;
        elspec.contact_diameter=0.8;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.8;
        elspec.tip_color=0.3;
        elspec.tip_length=2;
        elspec.contact_spacing=1.5;
        elspec.numel=18;
        elspec.tipiscontact=1;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)','K14 (R)','K15 (R)','K16 (R)','K17 (R)',...
            'K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (L)','K27 (L)','K28 (L)','K29 (L)','K30 (L)','K31 (L)','K32 (L)','K33 (L)','K34 (L)','K35 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech SD10R-SP05X Choi'
        elspec.matfname='adtech_sd10r_sp05x_choi';
        elspec.lead_diameter=1.1;
        elspec.lead_color=0.7;
        elspec.contact_length=2.4;
        elspec.contact_diameter=1.1;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.1;
        elspec.tip_color=0.7;
        elspec.tip_length=1;
        elspec.contact_spacing=2.4;
        elspec.numel=10;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech RD10R-SP03X'
        elspec.matfname='adtech_rd10r_sp03x';
        elspec.lead_diameter=0.86;
        elspec.lead_color=0.7;
        elspec.contact_length=2.29;
        elspec.contact_diameter=0.86;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.86;
        elspec.tip_color=0.7;
        elspec.tip_length=1;
        elspec.contact_spacing=0.71;
        elspec.numel=10;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech BF08R-SP05X'
        elspec.matfname='adtech_bf08r_sp05x';
        elspec.lead_diameter=1.28;
        elspec.lead_color=0.7;
        elspec.contact_length=1.57;
        elspec.contact_diameter=1.28;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.28;
        elspec.tip_color=0.7;
        elspec.tip_length=1;
        elspec.contact_spacing=3.43;
        elspec.numel=8;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech BF08R-SP21X'
        elspec.matfname='adtech_bf08r_sp21x';
        elspec.lead_diameter=1.3;
        elspec.lead_color=0.7;
        elspec.contact_length=1.57;
        elspec.contact_diameter=1.3;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.3;
        elspec.tip_color=0.7;
        elspec.tip_length=1;
        elspec.contact_spacing={1.43, 3.93};
        elspec.numel=8;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech BF08R-SP61X'
        elspec.matfname='adtech_bf08r_sp61x';
        elspec.lead_diameter=1.3;
        elspec.lead_color=0.7;
        elspec.contact_length=1.57;
        elspec.contact_diameter=1.3;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.3;
        elspec.tip_color=0.7;
        elspec.tip_length=1;
        elspec.contact_spacing={1.43, 4.43};
        elspec.numel=8;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'ELAINE Rat Electrode'
        elspec.matfname='elaine_rat_electrode';
        elspec.lead_diameter=225/1000;
        elspec.lead_color=0.7;
        elspec.contact_length=0;
        elspec.contact_diameter=225/1000;
        elspec.contact_color=0.3;
        elspec.tip_diameter=225/1000;
        elspec.tip_color=0.7;
        elspec.tip_length=112.5/1000;
        elspec.contact_spacing=0;
        elspec.numel=1;
        elspec.tipiscontact=1;
        elspec.contactnames={'K0 (R)','K0 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'NuMed Mini Lead'
        elspec.matfname='numed_minilead';
        elspec.lead_diameter=0.3125;
        elspec.lead_color=0.7;
        elspec.contact_length=0.5;
        elspec.contact_diameter=0.3125;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.3125;
        elspec.tip_color=0.7;
        elspec.tip_length=0.5;
        elspec.contact_spacing=0.5;
        elspec.numel=4;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K0 (L)','K1 (L)','K2 (L)','K3 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
end

if ~isfield(elspec,'eldist') && isa(elspec.contact_spacing, 'cell')
    elspec.eldist=elspec.contact_spacing{2}+elspec.contact_length;
elseif ~isfield(elspec,'eldist')
    elspec.eldist=elspec.contact_spacing+elspec.contact_length;
end

try
    options.elspec=elspec;
catch
    keyboard
end

varargout{1}=options;
