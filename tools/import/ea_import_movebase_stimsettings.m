function ea_import_movebase_stimsettings(hobj, evt, handles)

patient_dir = getappdata(handles.leadfigure,'uipatdir');

% sanity check if only one patient has been selected
if length(patient_dir) > 1
    disp('More than one patient selected, but this import function only works for one patient!')
    return;
end

if ~isempty(patient_dir)
    options = ea_getptopts(patient_dir{1});
else
    disp('No patient loaded. Please select a patient first.')
    return;
end

% get the .json file and load data
[json_file, json_path] = uigetfile('*.json');
json_data = loadjson(fullfile(json_path, json_file));

procedure_data = get_procedure_data(json_data);

% sanity check if electrodes have actually been found, otherwise exit
if isempty(procedure_data)
    disp('No electrodes found in the .json file provided!');
    return;
end

% get stimulation amplitude, frequencey and pulse width for both hemispheres
general_stim_settings = get_stim_amp_freq_pwidth(procedure_data);

% get polarity and percentages of contacts and casings
contact_stim_settings = get_contact_settings(procedure_data);

% sanity check if number of electrodes equals to the one of the selected electrode
if ~(options.elspec.numContacts == length(contact_stim_settings{1, 1}.contact_polarity))
    disp('Number of contacts of selected electrode and imported .json file do not match!')
    return;
end

% get name for stim
stim_name = inputdlg('Please specify a stimulation name');

S_init = ea_initializeS(stim_name, options);    % init S

S = construct_S(S_init, contact_stim_settings, general_stim_settings);  % fill S values\

if strcmp(general_stim_settings{1, 1}.stim_amp.unit, 'Volt')
    S=calculate_even_voltage_percentage(S);  % get correct percentage if voltage controlled
end

ea_savestimulation(S,options)   % create stim folder and stimparameters.mat

fprintf('Successfully created stimulation at %s', fullfile(patient_dir{1}, 'stimulations', ea_getspace, S.label))
end

%% utility functions
function output = get_procedure_data(json_input)

% find the entries that are labelled as procedures, as each procedure is an electrode in FHIR/dotbase
output = {};
for i = 1:length(json_input.entry)
    if strcmp(json_input.entry{1, i}.resource.resourceType, 'Procedure')
        output{end + 1} = json_input.entry{1, i};
    end
end
end

function output = get_stim_amp_freq_pwidth(procedure_data)

output = {};
for hemi = 1:2

    for extension_nr = 1:length(procedure_data{1, hemi}.resource.extension{1, 1}.extension)
        if strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.url, 'amplitude')
            output{hemi}.stim_amp.value = procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.valueQuantity.value;
            output{hemi}.stim_amp.unit = procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.valueQuantity.unit;
        elseif strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.url, 'frequency')
            output{hemi}.freq.value = procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.valueQuantity.value;
            output{hemi}.freq.unit = procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.valueQuantity.unit;
        elseif strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.url, 'pulse-width')
            output{hemi}.pwidth.value = procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.valueQuantity.value;
            output{hemi}.pwidth.unit = procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.valueQuantity.unit;

        end
    end
end
end

function output = get_contact_settings(procedure_data)

output = {};
for hemi = 1:2

    output{hemi}.contact_polarity = {};
    output{hemi}.case_polarity = {};

    output{hemi}.contact_percentage = {};
    output{hemi}.case_percentage = {};


    for extension_nr = 1:length(procedure_data{1, hemi}.resource.extension{1, 1}.extension)

        % if extension is a field (so it is not amplitude, frequency or pulse-width
        if isfield(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}, 'extension')
            % contact but not case
            if strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 1}.url, 'contact-label') && ...
                    ~strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 1}.valueCoding.display, 'case')

                if strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 3}.valueCode, 'n')
                    output{hemi}.contact_polarity{end + 1, 1} = 1;
                elseif strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 3}.valueCode, 'p')
                    output{hemi}.contact_polarity{end + 1, 1} = 2;
                else
                    output{hemi}.contact_polarity{end + 1, 1} = 0;
                end

                %  get percentages of contacts
                if isfield(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 2}.valueQuantity, 'value')
                    output{hemi}.contact_percentage{end + 1, 1} = procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 2}.valueQuantity.value;
                else
                    output{hemi}.contact_percentage{end + 1, 1} = 0;
                end

                % case
            elseif strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 1}.url, 'contact-label') && ...
                    strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 1}.valueCoding.display, 'case')

                if strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 3}.valueCode, 'n')
                    output{hemi}.case_polarity{end + 1, 1} = 1;
                elseif strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 3}.valueCode, 'p')
                    output{hemi}.case_polarity{end + 1, 1} = 2;
                else
                    output{hemi}.case_polarity{end + 1, 1} = 0;
                end

                % get percentages of case

                if isfield(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 2}.valueQuantity, 'value')
                    output{hemi}.case_percentage{end + 1, 1} = procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 2}.valueQuantity.value;
                else
                    output{hemi}.case_percentage{end + 1, 1} = 0;
                end

            end
        end
    end

end
end

function output = construct_S(S_init, contact_stim_settings, general_stim_settings)

output = S_init;

for hemi = 1:2

    % set amplitudes
    output.amplitude{1, hemi}(1) = general_stim_settings{1, hemi}.stim_amp.value;

    % 1 - right hemisphere, 2 - left hemisphere
    if hemi == 1
        prefix = 'Rs';
    else
        prefix = 'Ls';
    end

    % set amplitudes (again)
    output.([prefix, '1']).amp = general_stim_settings{1, hemi}.stim_amp.value;

    if strcmp(general_stim_settings{1, hemi}.stim_amp.unit, 'Volt')
        output.([prefix, '1']).va = 1;      % voltage controlled - 1, current controlled - 2
    else
        output.([prefix, '1']).va = 2;      % voltage controlled - 1, current controlled - 2
    end

    % set case
    if ~contact_stim_settings{1, hemi}.case_polarity{1} == 0
        output.([prefix, '1']).case.perc = 100;
    end
    output.([prefix, '1']).case.pol = contact_stim_settings{1, hemi}.case_polarity{1};

    % set contacts
    for contactNr = 1:length(contact_stim_settings{1, hemi}.contact_polarity)
        output.([prefix, num2str(1)]).(['k', num2str(contactNr)]).pol = contact_stim_settings{1, hemi}.contact_polarity{contactNr, 1};
        output.([prefix, num2str(1)]).(['k', num2str(contactNr)]).perc = contact_stim_settings{1, hemi}.contact_percentage{contactNr, 1};
    end
end

end

function S = calculate_even_voltage_percentage(S)

for hemi = 1:2

    % 1 - right hemisphere, 2 - left hemisphere
    if hemi == 1
        prefix = 'Rs';
    else
        prefix = 'Ls';
    end

    % find out how many contacts have positive and negative contacts ( 1- negative, 2 - positive)
    nr_pos_contacts = 0;
    nr_neg_contacts = 0;
    for contactNr = 1:S.numContacts
        if S.([prefix, '1']).(['k', num2str(contactNr)]).pol == 1
            nr_neg_contacts = nr_neg_contacts + 1;
        elseif S.([prefix, '1']).(['k', num2str(contactNr)]).pol == 2
            nr_pos_contacts = nr_pos_contacts + 1;
        end
    end

    % do the same for case
    if S.([prefix, '1']).case.pol == 1
        nr_neg_contacts = nr_neg_contacts + 1;
    elseif S.([prefix, '1']).case.pol == 2
        nr_pos_contacts = nr_pos_contacts + 1;
    end

    % small sanity check if no positive/negative contacts have been found
    if nr_neg_contacts == 0
        fprintf('No negative contacts found, please assign negative contacts manually after import.')
        nr_neg_contacts = 1;   % set number to 1 to avoid division by 0
    elseif nr_pos_contacts == 0
        fprintf('No positive contacts found, please assign negative contacts manually after import.')
        nr_pos_contacts = 1;   % set number to 1 to avoid division by 0
    end

    % now go through them again and assign the correct percentage
    for contactNr = 1:S.numContacts
        if S.([prefix, '1']).(['k', num2str(contactNr)]).pol == 1
            S.([prefix, '1']).(['k', num2str(contactNr)]).perc = 100 / nr_neg_contacts;
        elseif S.([prefix, '1']).(['k', num2str(contactNr)]).pol == 2
            S.([prefix, '1']).(['k', num2str(contactNr)]).perc = 100 / nr_pos_contacts;
        end
    end

    % do the same for case
    if S.([prefix, '1']).case.pol == 1
        S.([prefix, '1']).case.perc = 100 / nr_neg_contacts;
    elseif S.([prefix, '1']).case.pol == 2
        S.([prefix, '1']).case.perc = 100 / nr_pos_contacts;
    end

end
end
