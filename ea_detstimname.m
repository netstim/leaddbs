function [stimname,preexist]=ea_detstimname(options)
preexist = 0;

% check if previous stimulations have been stored
if ~isfield(options,'root') % called from lead group
    stimname = ['gs_',options.groupid];
    return
end

stimname = cell(0);
if exist([options.subj.stimDir, filesep, ea_nt(options)],'dir')
    stimdir = dir([options.subj.stimDir, filesep, ea_nt(options)]);
    [~, ind] = sort([stimdir(:).datenum], 'descend'); % show the latest modified first
    stimdir = stimdir(ind);
    stimname = {stimdir(cell2mat({stimdir.isdir})).name};
    stimname = stimname(cellfun(@(x) ~startsWith(x, {'.', 'gs_'}), stimname));
end

if ~isempty(stimname)
    preexist = 1;
end

if isempty(stimname) || (isfield(options, 'gen_newstim') && options.gen_newstim==1)
    stimname{end+1} = ea_getnewstimname;
end

% add commands
stimname{end+1} = ' => New stimulation';
stimname{end+1} = ' => Rename stimulation';
stimname{end+1} = ' => Delete stimulation';
stimname = stimname';


function stimname=ea_getnewstimname
try
    stimname=datestr(datevec(now), 'yyyymmddHHMMSS' );
catch
    import java.util.UUID;

    uid = char(UUID.randomUUID());
end
while 1
    stimc = inputdlg('Please enter a label for this stimulation','Stimulation Label',1,{stimname});
    if length(stimc{1})<3
        break
    else
        if startsWith(stimc{1},'gs_')
            msgbox('Please do not choose a stimulation label that starts with "gs_". These are reserved letters used in stimulations programmed inside Lead Group.')
        else
            break
        end
    end
end
stimname=stimc{1};
