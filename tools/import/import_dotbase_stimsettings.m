function import_dotbase_stimsettings(hobj, evt, handles)

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
contact_stim_settings = get_contact_settings(procedure_data, general_stim_settings{1, 1}.stim_amp.unit);

% sanity check if number of electrodes equals to the one of the selected electrode
if ~(options.elspec.numel == length(contact_stim_settings{1, 1}.contact_polarity))
    disp('Number of contacts of selected electrode and imported .json file do not match!')
    return;
end

% get name for stim
stim_name = inputdlg('Please specify a stimulation name');

S_init = ea_initializeS(stim_name, options);    % init S

S = construct_S(S_init, contact_stim_settings, general_stim_settings);  % fill S values\

if strcmp(general_stim_settings{1, 1}.stim_amp.unit, 'Volt')
    S=ea_redistribute_voltage(S, options);  % get correct percentage if voltage controlled
end

ea_savestimulation(S,options)   % create stim folder and stimparameters.mat

fprintf('Successfully created stimulation at %s', fullfile())
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
        contact_offset = -1;
    else
        prefix = 'Ls';
        contact_offset = 7;
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
        output.([prefix, num2str(1)]).(['k', num2str(contactNr + contact_offset)]).pol = contact_stim_settings{1, hemi}.contact_polarity{contactNr, 1};
        output.([prefix, num2str(1)]).(['k', num2str(contactNr + contact_offset)]).perc = contact_stim_settings{1, hemi}.contact_percentage{contactNr, 1};
    end
end

end

function S=ea_redistribute_voltage(S,changedobj)
Rconts={'k0','k1','k2','k3','k4','k5','k6','k7'};
Lconts={'k8','k9','k10','k11','k12','k13','k14','k15'};
LcontsCase=[Lconts,{'case'}];
RcontsCase=[Rconts,{'case'}];
if ischar(changedobj) % different polarity on the block
    switch changedobj
        case Rconts
            conts=Rconts;
            contsCase=RcontsCase;
            sidec='R';
            side=1;
        case Lconts
            conts=Lconts;
            contsCase=LcontsCase;
            sidec='L';
            side=2;
        case 'Rcase'
            conts=Rconts;
            changedobj='case';
            contsCase=RcontsCase;

            side=1;
            sidec='R';
        case 'Lcase'
            conts=Lconts;
            contsCase=LcontsCase;

            changedobj='case';
            side=2;
            sidec='L';
    end

    % check polarity of changed object:
    polchanged=eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol']);

    % check for monopolar models:
    if S.monopolarmodel % these allow only 1 active anode contact per model.
        for c=1:length(conts)
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.pol=0;']);
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.perc=0;']);
        end
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol=1;']);
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100;']);

        return

    else
        %         if S.([sidec,'s',num2str(S.active(side))]).va==2 % ampere only allows one anode and one cathode
        %             for c=1:length(contsCase)
        %
        %                 if S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol==polchanged % same polarity as changed object
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol=ea_swappol(polchanged);
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).perc=100;
        %                 else
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol=0;
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).perc=0;
        %                 end
        %             end
        %             S.([sidec,'s',num2str(S.active(side))]).(changedobj).pol=1;
        %             S.([sidec,'s',num2str(S.active(side))]).(changedobj).perc=100;
        %         end
    end

    if polchanged==0
        % set changed contacts percentage to zero:
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=0;']);
    else
        % determine how many other nodes with this polarity exist:
        divby=1;
        contacts={};
        for con=1:length(conts)
            if eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.pol==polchanged'])
                if ~strcmp(conts{con},changedobj)
                    %voltages{divby}=eval(['S.Rs',num2str(S.active(side)),'.',Rconts{con},'.perc']);
                    contacts{divby}=conts{con};
                    divby=divby+1;
                end
            end
        end

        if eval(['S.',sidec,'s',num2str(S.active(side)),'.case.pol==polchanged'])
            if ~strcmp(changedobj,'case')
                contacts{divby}='case';
                divby=divby+1;
            end
        end
        % add case to calculation.

        % set changed contacts percentage:
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100/divby;']);

        % reduce all other contacts percentages:

        try divby=divby/length(contacts); end
        for c=1:length(contacts)
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc=',...
                'S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc/divby;']);
        end
    end

    % now clean up mess from polarity that the contact used to have..

    polchanged=ea_polminus(polchanged);
    sumpercs=0;

    if polchanged % polarization has changed from negative to positive. clean up negatives. or changed from positive to off. clean up positives.
        contacts={};
        cnt=0;
        for con=1:length(conts)
            if eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.pol==polchanged'])
                if ~strcmp(conts{con},changedobj)
                    %voltages{divby}=eval(['S.Rs',num2str(S.active(side)),'.',Rconts{con},'.perc']);

                    cnt=cnt+1;
                    contacts{cnt}=conts{con};
                    sumpercs=sumpercs+eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.perc']);
                end
            end
        end
        % add case to calculation:
        if eval(['S.',sidec,'s',num2str(S.active(side)),'.case.pol==polchanged'])
            if ~strcmp(changedobj,'case')
                cnt=cnt+1;
                contacts{cnt}='case';
                sumpercs=sumpercs+eval(['S.',sidec,'s',num2str(S.active(side)),'.case.perc']);
            end
        end

        multby=(100/sumpercs);
        if cnt
            for c=1:length(contacts)
                eval(['S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc=',...
                    'S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc*multby;']);
            end
        end
    end

else % voltage percentage changed
    changedobj=get(changedobj,'Tag');
    changedobj=changedobj(1:end-1);

    switch changedobj
        case Rconts
            conts=Rconts;
            sidec='R';
            side=1;
        case Lconts
            conts=Lconts;
            sidec='L';
            side=2;
        case 'RC'
            conts=Rconts;
            changedobj='case';
            side=1;
            sidec='R';
        case 'LC'
            conts=Lconts;
            changedobj='case';
            side=2;
            sidec='L';
    end

    % check for monopolar models:
    if S.monopolarmodel % these allow only 1 active anode contact per model.
        for c=1:length(conts)
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.pol=0;']);
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.perc=0;']);
        end
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol=1;']);
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100;']);

        return
    end

    % check polarity of changed object:
    try
        polchanged=eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol']);
    catch
        keyboard
    end

    if polchanged==0 % set changed contacts polarity to negative
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol=1;']);
        polchanged=1;
    end

    % determine how many other nodes with this polarity exist:
    divby=1;
    contacts={};
    sumpercent=0;
    for con=1:length(conts)
        if eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.pol==polchanged'])
            if ~strcmp(conts{con},changedobj)
                sumpercent=sumpercent+eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.perc']);
                contacts{divby}=conts{con};
                divby=divby+1;
            end
        end
    end

    % add case to calculation.
    if eval(['S.',sidec,'s',num2str(S.active(side)),'.case.pol==polchanged'])
        if ~strcmp(changedobj,'case')
            contacts{divby}='case';
            divby=divby+1;
        end
    end

    if divby==1 % only one contact -> set to 100 percent.
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100;']);
    end

    % reduce all other contacts percentages:
    divby=sumpercent/(100-eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc']));

    for c=1:length(contacts)
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc=',...
            'S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc/divby;']);
    end

end