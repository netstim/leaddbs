function emagn = ea_get_ossdbs_StimSets_VTA(stim_vector, efield_contact_solutions, stim_i)
% Superpostion of Lead-DBS BIDS format E-field solutions for unit contact-ground solutions.
% Stores superimposed 4-D Efields (in V/mm)
% Also stores E-field magnitude separately (in V/m)
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    stim_vector                 % 1xN vector, stimulation current over contacts (in A!)
    efield_contact_solutions    % Nx1 cell, niftis of unit contact-ground solutions
    stim_i                      % index of the stimulation protocol
end

% initialize and scale for the first contact
e_field_superimposed = efield_contact_solutions{1,1};
% assume output name based on the unit solutions and stim_i
e_field_superimposed.fname = [efield_contact_solutions{1,1}.fname(1:end-6),'_protocol_',num2str(stim_i),'.nii'];
e_field_superimposed.img(:) = e_field_superimposed.img(:) * stim_vector(1,1);

% superimpose the rest contacts
for contact_i = 2:length(stim_vector)
    e_field_superimposed.img(:) = e_field_superimposed.img(:) + efield_contact_solutions{contact_i,1}.img(:) * stim_vector(1,contact_i);
end


% also save 3-D nifti with magnitude
ea_dispt('Creating nifti header for export...');
% create nifti
[~, ~, endian] = computer;
switch endian
    case 'L'
        endian = 0;
    case 'B'
        endian = 1;
end
[stim_folder,~,~] = fileparts(efield_contact_solutions{1,1}.fname);
emagn_img = sqrt(dot(e_field_superimposed.img,e_field_superimposed.img,4))*1000.0;
emagn.img = emagn_img;
emagn.mat = e_field_superimposed.mat;
emagn.dim = e_field_superimposed.dim;
emagn.dt  = [4, endian];
emagn.n=[1 1];
emagn.descrip='oss-dbs-v2 - StimSets Magn';
return