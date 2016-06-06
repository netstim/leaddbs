function stimname=ea_detstimname()
try
    stimname=datestr(datevec(now), 'yyyymmddHHMMSS' );
catch
    import java.util.UUID;
    
    uid = char(UUID.randomUUID());
end

stimc = inputdlg('Please enter a label for this stimulation','Stimulation Label',1,{stimname});
stimname=stimc{1};