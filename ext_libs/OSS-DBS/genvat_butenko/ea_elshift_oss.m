function settings = ea_elshift_oss(settings,side,first_active,last_active)
% By Butenko, konstantinmgtu@gmail.com
% OSS assumes straight electrodes building them from the first
% contact. So if active contacts are too distal, we need to
% imitate electrode shift.
% IMPORTANT: this trick only works if the active contacts are
% not far away from each other!

arguments
    settings            % parameters for OSS-DBS simulation
    side                {mustBeNumeric} % hemisphere index (0 - rh, 1 - lh)
    first_active        {mustBeNumeric} % index of the first active contact on the lead (count starts from 1)
    last_active         {mustBeNumeric} % index of the last active contact on the lead (count starts from 1)
end    

% we will not shift to the first contact, since it might be
% active tip. 
% we need to take into account that tip might have a different length
tip_center_to_C2 = norm(settings.contactLocation{side}(2,:) - settings.contactLocation{side}(1,:));
% for other contacts it is assumed the same
intercontact_spacing = norm(settings.contactLocation{side}(first_active,:) - settings.contactLocation{side}(first_active-1,:));
% we enforce straight trajectory between one before first
% active contact and the last active!
traj_unit_vector = (settings.contactLocation{side}(last_active,:) - settings.contactLocation{side}(first_active-1,:)) / norm(settings.contactLocation{side}(last_active,:) - settings.contactLocation{side}(first_active-1,:));

settings.Implantation_coordinate(side,:) = settings.contactLocation{side}(first_active,:) - tip_center_to_C2 * traj_unit_vector;
% robust implementation for tail
if strcmp(settings.Electrode_type,"Boston Scientific Vercise")
    settings.Second_coordinate(side,:) = settings.contactLocation{side}(first_active-1,:) + tip_center_to_C2 * traj_unit_vector + 6 * intercontact_spacing * traj_unit_vector;
else
    settings.Second_coordinate(side,:) = settings.contactLocation{side}(first_active-1,:) + tip_center_to_C2 * traj_unit_vector + 2 * intercontact_spacing * traj_unit_vector;
end

% we also need to adjust stimulation
settings.Phi_vector_old = settings.Phi_vector;
% per definition for shifted lead, 2nd contact is first active
settings.Phi_vector(side,:) = nan;
settings.Phi_vector(side,2:2+last_active-first_active) = settings.Phi_vector_old(side,first_active:last_active);
