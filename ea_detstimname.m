function stimname=ea_detstimname()
try
    stimname=datestr(datevec(now), 'yyyy-mm-dd HH:MM:SS' );
catch
    import java.util.UUID;
    
    uid = char(UUID.randomUUID());
end