function electrode=ea_elspec_stjude_directed_15(varargin)

% This function creates the electrode specification for the short SJM directed lead.
% In contrast to other electrode generation functions, it is based on an
% externally generated FEM-compatible model, stored in the
% .\StJude_Directed_05_Components\ subfolder.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

electrodeorder = [1 2 4 3 5 7 6 8 9]; % 211117 - small change so that order of the directional electrodes is clockwise seen from the tip

%% import insulations and contacts from subfolder
for k = 1:18
    filename = ['.\StJude_Directed_15_Components\Insulations\' 'ins' num2str(k) '.1'];
    [node,~,face]=readtetgen(filename);
    electrode.insulation(k).vertices = node;
    electrode.insulation(k).faces = face(:,1:3);
    clear face node filename
end

for k = 1:9
%     filename = ['.\Boston_Vercise_Directed_Components\Contacts\' 'con' num2str(k) '.1'];
    filename = ['.\StJude_Directed_15_Components\Contacts\' 'con' num2str(electrodeorder(k)) '.1'];
    [node,~,face]=readtetgen(filename);
    electrode.contacts(k).vertices = node;
    electrode.contacts(k).faces = face(:,1:3);
    clear face node filename
end


%% other specifications
electrode.electrode_model = 'St. Jude Directed 6173 (long)';
electrode.head_position = [0 0 0.75];
electrode.tail_position = [0 0 6.75];
electrode.x_position = [0.65 0 0.75];
electrode.y_position = [0 0.65 0.75];
electrode.numel = 8;
electrode.contact_color = 0.3;
electrode.lead_color = 0.7;

electrode.coords_mm(1,:)=[0 0 0.75];
electrode.coords_mm(2,:)=[0 0 3.75]+[-0.66,0,0];
electrode.coords_mm(3,:)=[0 0 3.75]+[0.33,0.66,0];
electrode.coords_mm(4,:)=[0 0 3.75]+[0.33,-0.66,0];
electrode.coords_mm(5,:)=[0 0 6.75]+[-0.66,0,0];
electrode.coords_mm(6,:)=[0 0 6.75]+[0.33,0.66,0];
electrode.coords_mm(7,:)=[0 0 6.75]+[0.33,-0.66,0];
electrode.coords_mm(8,:)=[0 0 9.75];


%% saving electrode struct
save('stjude_directed_15.mat','electrode');

%% create and save _vol file
filename = ['.\StJude_Directed_15_Components\final.1'];
[node,~,face]=readtetgen(filename);
save('stjude_directed_15_vol.mat','face','node')
clear node face
end