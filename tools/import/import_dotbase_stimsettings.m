function import_dotbase_stimsettings()

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

contact_stim_settings = get_contact_settings(procedure_data, general_stim_settings{1, 1}.stim_amp.unit);


disp('done')

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


function output = get_contact_settings(procedure_data, stim_unit)

output = {};
for hemi = 1:2

    output{hemi}.contact_polarity = {};
    output{hemi}.case_polarity = {};

    % if current controlled, add percentage
    if ~strcmp(stim_unit, 'Volt')
        output{hemi}.contact_percentage = {};
    end

    for extension_nr = 1:length(procedure_data{1, hemi}.resource.extension{1, 1}.extension)

        % if extension is a field (so it is not amplitude, frequency or pulse-width
        if isfield(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}, 'extension')
            % contact but not case
            if strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 1}.url, 'contact-label') && ~strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 1}.valueCoding.display, 'case')
                if strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 3}.valueCode, 'n')
                    output{hemi}.contact_polarity{end + 1, 1} = -1;
                elseif strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 3}.valueCode, 'p')
                    output{hemi}.contact_polarity{end + 1, 1} = 1;
                else
                    output{hemi}.contact_polarity{end + 1, 1} = 0;
                end
                % case
            elseif strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 1}.url, 'contact-label') && strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 1}.valueCoding.display, 'case')
                if strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 3}.valueCode, 'n')
                    output{hemi}.case_polarity{end + 1, 1} = -1;
                elseif strcmp(procedure_data{1, hemi}.resource.extension{1, 1}.extension{1, extension_nr}.extension{1, 3}.valueCode, 'p')
                    output{hemi}.case_polarity{end + 1, 1} = 1;
                else
                    output{hemi}.case_polarity{end + 1, 1} = 0;
                end
            end
        end
    end

end
end