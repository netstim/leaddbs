function stimname=ea_detstimname(options,handles)

% check if previous stimulations have been stored
if ~isfield(options,'root') % called from lead group
   stimname='gs'; 
   return
end
directory=[options.root,options.patientname,filesep];
stimname=cell(0);
if exist([directory,'stimulations'],'dir')
    sd=dir([directory,'stimulations']);
    for s=1:length(sd)
       if sd(s).isdir && ~strcmp(sd(s).name(1),'.')
           stimname{end+1}=sd(s).name;
       end
        
    end
    
end


if isempty(stimname) || (isfield(options,'gen_newstim') && options.gen_newstim==1)
    stimname{end+1}=ea_getnewstimname;
end

% add commands
stimname{end+1}=' => New stimulation';
stimname{end+1}=' => Rename stimulation';
stimname{end+1}=' => Delete stimulation';
stimname=stimname';



function stimname=ea_getnewstimname
try
    stimname=datestr(datevec(now), 'yyyymmddHHMMSS' );
catch
    import java.util.UUID;
    
    uid = char(UUID.randomUUID());
end

stimc = inputdlg('Please enter a label for this stimulation','Stimulation Label',1,{stimname});
stimname=stimc{1};