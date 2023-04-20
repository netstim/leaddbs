function [contactnames,directional]=ea_getelcontactnames(model,side)
% generic / lead-dbs:
switch side
    case 1 % RH
        contactnames={'k0','k1','k2','k3','k4','k5','k6','k7'};
        directional=[0,0,0,0,0,0,0,0];
    case 2 % LH
        contactnames={'k7','k8','k9','k10','k11','k12','k13','k14'};
        directional=[0,0,0,0,0,0,0,0];
end



switch model
    case {'Medtronic 3389','Medtronic 3387','Medtronic 3391'} % classics
        switch side
            case 1 % RH
                contactnames={'k0','k1','k2','k3'};
                directional=[0,0,0,0];
            case 2 % LH
                contactnames={'k7','k8','k9','k10'};
                directional=[0,0,0,0];
        end
    case {'Medtronic B33005','Medtronic B33015'} % sensight
        switch side
            case 1 % RH
                contactnames={'R0','R1A','R1B','R1C','R2A','R2B','R2C','R3'};
                directional=[0,1,1,1,1,1,1,0];
            case 2 % LH
                contactnames={'L0','L1A','L1B','L1C','L2A','L2B','L2C','L3'};
                directional=[0,1,1,1,1,1,1,0];
        end
    case 'Boston Scientific Vercise'
        switch side
            case 1 % RH
                contactnames={'R1','R2','R3','R4','R5','R6','R7','R8'};
                directional=[0,0,0,0,0,0,0,0];
            case 2 % LH
                contactnames={'L1','L2','L3','L4','L5','L6','L7','L8'};
                directional=[0,0,0,0,0,0,0,0];
        end
    case 'Boston Scientific Vercise Directed'
        switch side
            case 1 % RH
                contactnames={'R1','R2','R3','R4','R5','R6','R7','R8'};
                directional=[0,1,1,1,1,1,1,0];
            case 2 % LH
                contactnames={'L1','L2','L3','L4','L5','L6','L7','L8'};
                directional=[0,1,1,1,1,1,1,0];
        end
    case 'Boston Scientific Vercise Cartesia HX'
        switch side
            case 1 % RH
                contactnames={'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16'};
                directional=[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0];
            case 2 % LH
                contactnames={'L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12','L13','L14','L15','L16'};
                directional=[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0];
        end
    case 'Boston Scientific Vercise Cartesia X'
        switch side
            case 1 % RH
                contactnames={'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15','R16'};
                directional=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0];
            case 2 % LH
                contactnames={'L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12','L13','L14','L15','L16'};
                directional=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0];
        end
    case {'Abbott ActiveTip (6146-6149)','Abbott ActiveTip (6142-6145)'}
        switch side
            case 1 % RH
                contactnames={'R1','R2','R3','R4'};
                directional=[0,0,0,0];
            case 2 % LH
                contactnames={'L1','L2','L3','L4'};
                directional=[0,0,0,0];
        end
    case {'Abbott Directed 6172 (short)','Abbott Directed 6173 (long)'}
        switch side
            case 1 % RH
                contactnames={'R1','R2A','R2B','R2C','R3A','R3B','R3C','R4'};
                directional=[0,1,1,1,1,1,1,0];
            case 2 % LH
                contactnames={'L1','L2A','L2B','L2C','L3A','L3B','L3C','L4'};
                directional=[0,1,1,1,1,1,1,0];
        end
    case 'PINS Medical L301'
    case 'PINS Medical L302'
    case 'PINS Medical L303'
    case 'SceneRay SR1200'
    case 'SceneRay SR1210'
    case 'SceneRay SR1211'
    case 'SceneRay SR1242'
    case 'SDE-08 S8 Legacy'
    case 'SDE-08 S10 Legacy'
    case 'SDE-08 S12 Legacy'
    case 'SDE-08 S16 Legacy'
    case 'SDE-08 S8'
    case 'SDE-08 S10'
    case 'SDE-08 S12'
    case 'SDE-08 S14'
    case 'SDE-08 S16'
    case 'PMT 2102-16-092'
    case 'PMT 2102-16-093'
    case 'PMT 2102-16-131'
    case 'PMT 2102-16-142'
    case '2069-EPC-05C-35'
    case '2069-EPC-15C-35'
    case 'NeuroPace DL-344-3.5'
    case 'NeuroPace DL-344-10'
    case 'DIXI D08-05AM'
    case 'DIXI D08-08AM'
    case 'DIXI D08-10AM'
    case 'DIXI D08-12AM'
    case 'DIXI D08-15AM'
    case 'DIXI D08-18AM'
    case 'AdTech SD10R-SP05X Choi'
    case 'AdTech RD10R-SP03X'
    case 'AdTech BF08R-SP05X'
    case 'AdTech BF08R-SP21X'
    case 'AdTech BF08R-SP61X'
    case 'AdTech SD08R-SP05X'
    case 'AdTech SD14R-SP05X'
    case 'ELAINE Rat Electrode'
    case 'FHC WU Rat Electrode'
    case 'NuMed Mini Lead'
    case 'Aleva directSTIM Directed'
    case 'SmartFlow Cannula NGS-NC-06'
end

