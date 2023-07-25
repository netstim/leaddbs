function varargout=ea_resolve_elspec(varargin)
% This simple function outputs a cellarray of available electrode specs if
% nargin==0 and exports the current electrode specification if varargin{1}
% is an options struct with options.elmodel defined as a valid name of electrode string.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ~nargin
    varargout{1}={'Medtronic 3389', 'Medtronic 3387', 'Medtronic 3391', 'Medtronic B33005', 'Medtronic B33015', ...
        'Boston Scientific Vercise', 'Boston Scientific Vercise Directed', ...
        'Boston Scientific Vercise Cartesia HX', 'Boston Scientific Vercise Cartesia X', ...
        'Abbott ActiveTip (6146-6149)','Abbott ActiveTip (6142-6145)', ...
        'Abbott Directed 6172 (short)','Abbott Directed 6173 (long)', ...
        'PINS Medical L301', 'PINS Medical L302', 'PINS Medical L303', .....
        'SceneRay SR1200', 'SceneRay SR1210', 'SceneRay SR1211', 'SceneRay SR1242', ...
        'SDE-08 S8 Legacy', 'SDE-08 S10 Legacy', 'SDE-08 S12 Legacy', 'SDE-08 S16 Legacy', ...
        'SDE-08 S8', 'SDE-08 S10', 'SDE-08 S12', 'SDE-08 S14', 'SDE-08 S16', ...
        'PMT 2102-16-092', 'PMT 2102-16-093', 'PMT 2102-16-131', 'PMT 2102-16-142', ...
        '2069-EPC-05C-35', '2069-EPC-15C-35', 'NeuroPace DL-344-3.5', 'NeuroPace DL-344-10', ...
        'DIXI D08-05AM', 'DIXI D08-08AM', 'DIXI D08-10AM', 'DIXI D08-12AM', 'DIXI D08-15AM', 'DIXI D08-18AM', ...
        'AdTech BF08R-SP05X', 'AdTech BF08R-SP21X', 'AdTech BF08R-SP61X', 'AdTech BF09R-SP61X-0BB', ...
        'AdTech RD06R-SP05X', 'AdTech RD08R-SP05X', 'AdTech RD10R-SP03X', 'AdTech RD10R-SP05X', 'AdTech RD10R-SP06X', 'AdTech RD10R-SP07X', 'AdTech RD10R-SP08X', ...
        'AdTech SD06R-SP26X', 'AdTech SD08R-SP05X', 'AdTech SD10R-SP05X', 'AdTech SD10R-SP05X Choi', 'AdTech SD14R-SP05X', ...
        'ELAINE Rat Electrode', 'FHC WU Rat Electrode', 'NuMed Mini Lead', ...
        'Aleva directSTIM Directed', ...
        'SmartFlow Cannula NGS-NC-06'}';
    varargout{2}={'medtronic_3389', 'medtronic_3387', 'medtronic_3391', 'medtronic_b33005', 'medtronic_b33015', ...
        'boston_vercise', 'boston_vercise_directed', ...
        'boston_vercise_cartesia_hx', 'boston_vercise_cartesia_x', ...
        'abbott_activetip_2mm','abbott_activetip_3mm', ...
        'abbott_directed_05','abbott_directed_15', ...
        'pins_l301', 'pins_l302', 'pins_l303', ...
        'sceneray_sr1200', 'sceneray_sr1210', 'sceneray_sr1211', 'sceneray_sr1242', ...
        'sde_08_s8_legacy', 'sde_08_s10_legacy', 'sde_08_s12_legacy', 'sde_08_s16_legacy',...
        'sde_08_s8', 'sde_08_s10', 'sde_08_s12', 'sde_08_s14', 'sde_08_s16', ...
        'pmt_2102_16_092', 'pmt_2102_16_093', 'pmt_2102_16_131', 'pmt_2102_16_142', ...
        'epc_05c', 'epc_15c', 'neuropace_dl_344_35', 'neuropace_dl_344_10', ...
        'dixi_d08_05am', 'dixi_d08_08am', 'dixi_d08_10am', 'dixi_d08_12am', 'dixi_d08_15am', 'dixi_d08_18am', ...
        'adtech_bf08r_sp05x', 'adtech_bf08r_sp21x', 'adtech_bf08r_sp61x', 'adtech_bf09r_sp61x_0bb', ...
        'adtech_rd06r_sp05x', 'adtech_rd08r_sp05x', 'adtech_rd10r_sp03x', 'adtech_rd10r_sp05x', 'adtech_rd10r_sp06x', 'adtech_rd10r_sp07x', 'adtech_rd10r_sp08x', ...
        'adtech_sd08r_sp26x', 'adtech_sd08r_sp05x',  'adtech_sd10r_sp05x', 'adtech_sd10r_sp05x_choi', 'adtech_sd14r_sp05x', ...
        'elaine_rat_electrode', 'fhc_wu_rat_electrode', 'numed_minilead', ...
        'aleva_directstim_directed', ...
        'smartflow_ngs-nc-06'}';
    return
else
    options=varargin{1};
end

if ~isfield(options, 'elmodel')
    try
        load([options.root,options.patientname,filesep,'reconstruction',filesep,options.patientname,'_desc-reconstruction.mat'],'reco');
        elmodel = ea_get_first_notempty_elmodel(reco.props);
    catch
        %no model was found
        warning('No electrode model specified. Using Medtronic 3389.');
        elmodel = 'Medtronic 3389';
    end
else
    elmodel = options.elmodel;
end

elspec.elmodel = elmodel;

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
    case 'Medtronic B33005'
        elspec.matfname='medtronic_b33005';
        elspec.lead_diameter=1.3;
        elspec.lead_color=0.7;
        elspec.contact_length=1.5;
        elspec.contact_diameter=1.3;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.3;
        elspec.tip_color=0.7;
        elspec.tip_length=0.9;
        elspec.contact_spacing=0.5;
        elspec.numel=8;
        elspec.tipiscontact=0;
        elspec.markerpos = 14.4;
        elspec.markerlen = 2.3;
        elspec.contactnames={'K0 (R)','K1A (R)','K1B (R)','K1C (R)','K2A (R)','K2B (R)','K2C (R)','K3 (R)','K0 (L)','K1A (L)','K1B (L)','K1C (L)','K2A (L)','K2B (L)','K2C (L)','K3 (L)'};
        elspec.isdirected=1;
        elspec.etagenames{1}={'K0 (R)','K1 (R)','K2 (R)','K3 (R)'};
        elspec.etagenames{2}={'K0 (L)','K1 (L)','K2 (L)','K3 (L)'};
        elspec.etageidx={1,2:4,5:7,8};
        elspec.forstimulation=1;
    case 'Medtronic B33015'
        elspec.matfname='medtronic_b33015';
        elspec.lead_diameter=1.3;
        elspec.lead_color=0.7;
        elspec.contact_length=1.5;
        elspec.contact_diameter=1.3;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.3;
        elspec.tip_color=0.7;
        elspec.tip_length=0.9;
        elspec.contact_spacing=1.5;
        elspec.numel=8;
        elspec.tipiscontact=0;
        elspec.markerpos = 17.4;
        elspec.markerlen = 2.3;
        elspec.contactnames={'K0 (R)','K1A (R)','K1B (R)','K1C (R)','K2A (R)','K2B (R)','K2C (R)','K3 (R)','K0 (L)','K1A (L)','K1B (L)','K1C (L)','K2A (L)','K2B (L)','K2C (L)','K3 (L)'};
        elspec.isdirected=1;
        elspec.etagenames{1}={'K0 (R)','K1 (R)','K2 (R)','K3 (R)'};
        elspec.etagenames{2}={'K0 (L)','K1 (L)','K2 (L)','K3 (L)'};
        elspec.etageidx={1,2:4,5:7,8};
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
    case 'Boston Scientific Vercise Cartesia HX'
        elspec.matfname='boston_vercise_cartesia_hx';
        elspec.lead_diameter=1.3;
        elspec.lead_color=0.7;
        elspec.contact_length=1.5;
        elspec.contact_diameter=1.3;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.3;
        elspec.tip_color=0.3;
        elspec.tip_length=1.1;
        elspec.contact_spacing=0.5;
        elspec.numel=16;
        elspec.tipiscontact=0;
        elspec.markerpos = 20;
        elspec.markerlen = 3;
        elspec.contactnames={'K17 (R)','K18 (R)','K19 (R)','K20 (R)','K21 (R)','K22 (R)','K23 (R)','K24 (R)','K25 (R)','K26 (R)','K27 (R)','K28 (R)','K29 (R)','K30 (R)','K31 (R)','K32 (R)',...
            'K1 (L)','K2 (L)','K3 (L)','K4 (L)','K5 (L)','K6 (L)','K7 (L)','K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)'};
        elspec.isdirected=1;
        elspec.etagenames{1}={'K17-19 (R)','K20-22 (R)','K23-25 (R)','K26-28 (R)','K29 (R)','K30 (R)','K31 (R)','K32 (R)'};
        elspec.etagenames{2}={'K1-3 (L)','K4-6 (L)','K7-9 (L)','K10-12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)'};
        elspec.etageidx={1:3,4:6,7:9,10:12,13,14,15,16};
        elspec.forstimulation=1;
    case 'Boston Scientific Vercise Cartesia X'
        elspec.matfname='boston_vercise_cartesia_x';
        elspec.lead_diameter=1.3;
        elspec.lead_color=0.7;
        elspec.contact_length=1.5;
        elspec.contact_diameter=1.3;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.3;
        elspec.tip_color=0.3;
        elspec.tip_length=1.1;
        elspec.contact_spacing=0.5;
        elspec.numel=16;
        elspec.tipiscontact=0;
        elspec.markerpos = 16;
        elspec.markerlen = 3;
        elspec.contactnames={'K17 (R)','K18 (R)','K19 (R)','K20 (R)','K21 (R)','K22 (R)','K23 (R)','K24 (R)','K25 (R)','K26 (R)','K27 (R)','K28 (R)','K29 (R)','K30 (R)','K31 (R)','K32 (R)',...
            'K1 (L)','K2 (L)','K3 (L)','K4 (L)','K5 (L)','K6 (L)','K7 (L)','K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)'};
        elspec.isdirected=1;
        elspec.etagenames{1}={'K17-19 (R)','K20-22 (R)','K23-25 (R)','K26-28 (R)','K29-31 (R)','K32 (R)'};
        elspec.etagenames{2}={'K1-3 (L)','K4-6 (L)','K7-9 (L)','K10-12 (L)','K13-15 (L)','K16 (L)'};
        elspec.etageidx={1:3,4:6,7:9,10:12,13:15,16};
        elspec.forstimulation=1;
    case 'Abbott ActiveTip (6146-6149)'
        elspec.matfname='abbott_activetip_2mm';
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
    case 'Abbott ActiveTip (6142-6145)'
        elspec.matfname='abbott_activetip_3mm';
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
    case 'Abbott Directed 6172 (short)'
        elspec.matfname='abbott_directed_05';
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
    case 'Abbott Directed 6173 (long)'
        elspec.matfname='abbott_directed_15';
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
    case 'SceneRay SR1200'
        elspec.matfname='sceneray_sr1200';
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
        elspec.contactnames={'K4 (R)','K5 (R)','K6 (R)','K7 (R)','K0 (L)','K1 (L)','K2 (L)','K3 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'SceneRay SR1210'
        elspec.matfname='sceneray_sr1210';
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
        elspec.contactnames={'K4 (R)','K5 (R)','K6 (R)','K7 (R)','K0 (L)','K1 (L)','K2 (L)','K3 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'SceneRay SR1211'
        elspec.matfname='sceneray_sr1211';
        elspec.lead_diameter=1.27;
        elspec.lead_color=0.7;
        elspec.contact_length=1.5;
        elspec.contact_diameter=1.27;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.27;
        elspec.tip_color=0.7;
        elspec.tip_length=1.5;
        elspec.contact_spacing=1.0;
        elspec.numel=4;
        elspec.tipiscontact=0;
        elspec.contactnames={'K4 (R)','K5 (R)','K6 (R)','K7 (R)','K0 (L)','K1 (L)','K2 (L)','K3 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'SceneRay SR1242'
        elspec.matfname='sceneray_sr1242';
        elspec.lead_diameter=1.27;
        elspec.lead_color=0.7;
        elspec.contact_length=3.0;
        elspec.contact_diameter=1.27;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.27;
        elspec.tip_color=0.7;
        elspec.tip_length=1.5;
        elspec.contact_spacing=[2.0,4.0,4.0];
        elspec.numel=4;
        elspec.tipiscontact=0;
        elspec.contactnames={'K4 (R)','K5 (R)','K6 (R)','K7 (R)','K0 (L)','K1 (L)','K2 (L)','K3 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
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
        elspec.contact_spacing=4;
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
    case 'PMT 2102-16-131'
        elspec.matfname='pmt_2102_16_131';
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
        elspec.forstimulation=1;
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
        elspec.isdirected=0;
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
        elspec.lead_diameter=1.28;
        elspec.lead_color=0.7;
        elspec.contact_length=1.57;
        elspec.contact_diameter=1.28;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.28;
        elspec.tip_color=0.7;
        elspec.tip_length=1;
        elspec.contact_spacing=[1.43, 3.93, 3.93, 3.93, 3.93, 3.93, 3.93];
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
        elspec.lead_diameter=1.28;
        elspec.lead_color=0.7;
        elspec.contact_length=1.57;
        elspec.contact_diameter=1.28;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.28;
        elspec.tip_color=0.7;
        elspec.tip_length=1;
        elspec.contact_spacing=[1.43, 4.43, 4.43, 4.43, 4.43, 4.43, 4.43];
        elspec.numel=8;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech BF09R-SP61X-0BB'
        elspec.matfname='adtech_bf09r_sp61x_0bb';
        elspec.lead_diameter=1.28;
        elspec.lead_color=0.7;
        elspec.contact_length=1.57;
        elspec.contact_diameter=1.28;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.28;
        elspec.tip_color=0.7;
        elspec.tip_length=1;
        elspec.contact_spacing=[1.43, 4.43, 4.43, 4.43, 4.43, 4.43, 4.43, 4.43];
        elspec.numel=9;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)',...
            'K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech RD06R-SP05X'
        elspec.matfname='adtech_rd06r_sp05x';
        elspec.lead_diameter=0.86;
        elspec.lead_color=0.7;
        elspec.contact_length=2.29;
        elspec.contact_diameter=0.86;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.86;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=2.71;
        elspec.numel=6;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)',...
            'K6 (L)','K7 (L)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech RD08R-SP05X'
        elspec.matfname='adtech_rd08r_sp05x';
        elspec.lead_diameter=0.86;
        elspec.lead_color=0.7;
        elspec.contact_length=2.29;
        elspec.contact_diameter=0.86;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.86;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=2.71;
        elspec.numel=8;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)'};
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
        elspec.tip_length=2;
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
    case 'AdTech RD10R-SP05X'
        elspec.matfname='adtech_rd10r_sp05x';
        elspec.lead_diameter=0.86;
        elspec.lead_color=0.7;
        elspec.contact_length=2.29;
        elspec.contact_diameter=0.86;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.86;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=2.71;
        elspec.numel=10;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech RD10R-SP06X'
        elspec.matfname='adtech_rd10r_sp06x';
        elspec.lead_diameter=0.86;
        elspec.lead_color=0.7;
        elspec.contact_length=2.29;
        elspec.contact_diameter=0.86;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.86;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=3.71;
        elspec.numel=10;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech RD10R-SP07X'
        elspec.matfname='adtech_rd10r_sp07x';
        elspec.lead_diameter=0.86;
        elspec.lead_color=0.7;
        elspec.contact_length=2.29;
        elspec.contact_diameter=0.86;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.86;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=4.71;
        elspec.numel=10;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech RD10R-SP08X'
        elspec.matfname='adtech_rd10r_sp08x';
        elspec.lead_diameter=0.86;
        elspec.lead_color=0.7;
        elspec.contact_length=2.29;
        elspec.contact_diameter=0.86;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.86;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=5.71;
        elspec.numel=10;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech SD06R-SP26X'
        elspec.matfname='adtech_sd06r_sp26x';
        elspec.lead_diameter=1.12;
        elspec.lead_color=0.7;
        elspec.contact_length=2.41;
        elspec.contact_diameter=1.12;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.12;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=[2.59, 2.59, 7.59, 7.59, 7.59];
        elspec.numel=6;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)', ...
            'K6 (L)','K7 (L)','K8 (L)','K9 (L)','K10 (L)','K11 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech SD08R-SP05X'
        elspec.matfname='adtech_sd08r_sp05x';
        elspec.lead_diameter=1.12;
        elspec.lead_color=0.7;
        elspec.contact_length=2.41;
        elspec.contact_diameter=1.12;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.12;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=2.59;
        elspec.numel=8;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)',...
            'K8 (L)','K9 (L)','K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=1;
    case 'AdTech SD10R-SP05X'
        elspec.matfname='adtech_sd10r_sp05x';
        elspec.lead_diameter=1.12;
        elspec.lead_color=0.7;
        elspec.contact_length=2.41;
        elspec.contact_diameter=1.12;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.12;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=2.59;
        elspec.numel=10;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)',...
            'K10 (L)','K11 (L)','K12 (L)','K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)'};
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
    case 'AdTech SD14R-SP05X'
        elspec.matfname='adtech_sd14r_sp05x';
        elspec.lead_diameter=1.12;
        elspec.lead_color=0.7;
        elspec.contact_length=2.41;
        elspec.contact_diameter=1.12;
        elspec.contact_color=0.3;
        elspec.tip_diameter=1.12;
        elspec.tip_color=0.7;
        elspec.tip_length=2;
        elspec.contact_spacing=2.59;
        elspec.numel=14;
        elspec.tipiscontact=0;
        elspec.contactnames={'K0 (R)','K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)','K13 (R)',...
            'K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)','K25 (L)','K26 (R)','K27 (R)'};
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
    case 'FHC WU Rat Electrode'
        elspec.matfname='fhc_wu_rat_electrode';
        elspec.lead_diameter=125/1000;
        elspec.lead_color=0.7;
        elspec.contact_length=0;
        elspec.contact_diameter=125/1000;
        elspec.contact_color=0.3;
        elspec.tip_diameter=125/1000;
        elspec.tip_color=0.7;
        elspec.tip_length=250/1000;
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
    case 'Aleva directSTIM Directed'
        elspec.matfname='aleva_directstim_directed';
        elspec.lead_diameter=1.35;
        elspec.lead_color=0.7;
        elspec.contact_length=1.5;
        elspec.contact_diameter=1.35;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.5;
        elspec.tip_color=0.7;
        elspec.tip_length=1.1;
        elspec.contact_spacing=0.5;
        elspec.markerpos = 18; %TOCHECK
        elspec.numel=12;
        elspec.tipiscontact=0;
        elspec.contactnames={'K1 (R)','K2 (R)','K3 (R)','K4 (R)','K5 (R)','K6 (R)','K7 (R)','K8 (R)','K9 (R)','K10 (R)','K11 (R)','K12 (R)',...
            'K13 (L)','K14 (L)','K15 (L)','K16 (L)','K17 (L)','K18 (L)','K19 (L)','K20 (L)','K21 (L)','K22 (L)','K23 (L)','K24 (L)'};        
        elspec.isdirected=1;
        elspec.etagenames{1}={'K1-3 (R)','K4-6 (R)','K7-9 (R)','K10-12 (R)'};
        elspec.etagenames{2}={'K13-15 (L)','K16-18 (L)','K19-21 (L)','K22-24 (L)'};
        elspec.etageidx={1:3,4:6,7:9,10:12};
        elspec.forstimulation=1;
    case 'SmartFlow Cannula NGS-NC-06'
        elspec.matfname='smartflow_ngs-nc-06';
        elspec.lead_diameter=1.65;
        elspec.lead_color=0.7;
        elspec.contact_length=10;
        elspec.contact_diameter=0.27;
        elspec.contact_color=0.3;
        elspec.tip_diameter=0.2;
        elspec.tip_color=0.7;
        elspec.tip_length=3;
        elspec.contact_spacing=0;
        elspec.numel=2;
        elspec.tipiscontact=1;
        elspec.contactnames={'K0 (R)','K1 (R)','K3 (L)','K4 (L)'};
        elspec.isdirected=0;
        elspec.etagenames{1}=elspec.contactnames(1:length(elspec.contactnames)/2);
        elspec.etagenames{2}=elspec.contactnames((length(elspec.contactnames)/2)+1:end);
        elspec.etageidx=num2cell(1:elspec.numel);
        elspec.forstimulation=0;
end

if ~isfield(elspec,'eldist') && numel(elspec.contact_spacing)>1
    elspec.eldist=elspec.contact_spacing(2)+elspec.contact_length;
elseif ~isfield(elspec,'eldist')
    elspec.eldist=elspec.contact_spacing+elspec.contact_length;
end

try
    options.elspec=elspec;
catch
    keyboard
end

varargout{1}=options;
