%% Read Medtronic Framelink XML Configuration File
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function config = readMedtronicXMLTrajectory(filepath)
xml = xmlread(filepath);
config = struct();

%% Check Version
configVersion = char(xml.getElementsByTagName('dtd-version').item(0).getFirstChild.getData);
if(~(strcmp(configVersion, '1.0') || strcmp(configVersion, '2.0')))
    warning(['Only version 1.0 and 2.0 of the loaded Medtronic File is tested,  but the file you specified is Version ' configVersion]);
end

if(strcmp(configVersion, '2.0'))
    acPc = xml.getElementsByTagName('reformatSettings_v2').item(0);   
    surgicalPlans = xml.getElementsByTagName('surgicalPlan_v2');
    
    % Some Medtronic Files are buggy and use the old surgicalPlan tag
    % despite reporting v2.. try to fix this
    if(surgicalPlans.getLength() == 0)
        warning('.xml File reports Version 2.0 but contains no surgicalPlan_v2 tags. Either your file does not contain any plans or it is indicating the wrong version number. Trying to read old v1.0 tags...');
        surgicalPlans = xml.getElementsByTagName('surgicalPlan');
    end
    %% Read the reformat settings specific to Format V2.0
    m1Element = acPc.getElementsByTagName('midline1').item(0);
    config.M1 = get3DFloatCoordinate(m1Element);
    
    m2Element = acPc.getElementsByTagName('midline2').item(0);
    config.M2 = get3DFloatCoordinate(m2Element);
    
    m3Element = acPc.getElementsByTagName('midline3').item(0);
    config.M3 = get3DFloatCoordinate(m3Element);
    
elseif(strcmp(configVersion, '1.0'))
    acPc = xml.getElementsByTagName('ACPC').item(0);
    surgicalPlans = xml.getElementsByTagName('surgicalPlan');
    
    %% Read the reformat settings specific to Format V1.0
    m1Element = acPc.getElementsByTagName('midline').item(0);
    config.M1 = get3DFloatCoordinate(m1Element);
    
    m2Element = acPc.getElementsByTagName('midline').item(1);
    config.M2 = get3DFloatCoordinate(m2Element);
    
    m3Element = acPc.getElementsByTagName('midline').item(2);
    config.M3 = get3DFloatCoordinate(m3Element);
    
    frameRods = xml.getElementsByTagName('frameRods');
    frameRod = frameRods.item(0);
    
    for i=0:8
        rod = frameRod.getElementsByTagName('rod').item(i);
        config.rods(i+1).coord = get3DFloatCoordinate(rod);
    end
end

%% Read the reformat settings (== AC/PC/Midline defintions)
acElement = acPc.getElementsByTagName('AC').item(0);
config.AC = get3DFloatCoordinate(acElement);

pcElement = acPc.getElementsByTagName('PC').item(0);
config.PC = get3DFloatCoordinate(pcElement);


%% Read the surgical plans (== trajectories)
% all elementS (lists) that have a guaranted cardinality of 1 can be used with
% item(0) to give the respective element
noTraject = surgicalPlans.getLength();

for i = 0:noTraject-1;
    surgicalPlan = surgicalPlans.item(i);
    planNameElement = surgicalPlan.getElementsByTagName('name').item(0);
    config.trajects(i+1).name = char(planNameElement.getFirstChild.getData);
    
    targetElement = surgicalPlan.getElementsByTagName('target').item(0);
    config.trajects(i+1).target = get3DFloatCoordinate(targetElement);
    
    entryElement = surgicalPlan.getElementsByTagName('entry').item(0);
    config.trajects(i+1).entry  = get3DFloatCoordinate(entryElement);
end

%% Private Helper Functions
    function float3 = get3DFloatCoordinate(point3dActive)
        x = getFloatCoordinate(point3dActive, 'x');
        y = getFloatCoordinate(point3dActive, 'y');
        z = getFloatCoordinate(point3dActive, 'z');
        float3 = [x;y;z];
    end
    function float = getFloatCoordinate(point3dActive, strName)
        targetPoint3DElement = point3dActive.getElementsByTagName('point3dActive').item(0);
        coordElement = targetPoint3DElement.getElementsByTagName(strName).item(0);
        coordFloat = coordElement.getElementsByTagName('float').item(0);
        float = str2double(coordFloat.getFirstChild.getData);
    end
end