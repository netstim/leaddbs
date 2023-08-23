function ea_checkleaddirs
prefs=ea_prefs;

if ~isfolder(prefs.lc.datadir)
    if ~(isdeployed && ~isfield(prefs,'firstrun')) % don't enter first time in deployed app
        try
            ea_mkdir(prefs.lc.datadir);
        catch
            ea_cprintf('CmdWinWarnings', 'Failed to create connectome directory.\nPlease set prefs.lc.datadir to a valid directory in preferences.\n');
        end
    end
end

if ~isempty(prefs.ltx.pdfconverter)
   if ~isfile(prefs.ltx.pdfconverter)
        ea_cprintf('CmdWinWarnings', 'LaTeX PDF converter not set properly.\nPlease edit prefs.ltx.pdfconverter in preferences.\n'); 
   end
end

if ~isempty(prefs.ls.dir)
    if ~isfolder(prefs.ls.dir)
        try
            mkdir(prefs.ls.dir);
        catch
            ea_cprintf('CmdWinWarnings', 'Lead DBS Server directory could not be established.\nPlease set prefs.ls.dir folder to a valid directory in preferences.\n');
        end
        
    end
end

if ~isempty(prefs.ixi.dir)
    if ~isfolder(prefs.ixi.dir)
        try
            mkdir(prefs.ixi.dir);
        catch
            ea_cprintf('CmdWinWarnings', 'IXI data directory could not be established.\nPlease set prefs.ixi.dir folder to a valid directory in preferences.\n');
        end
    end
end
