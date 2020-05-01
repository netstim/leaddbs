function ea_checkleaddirs
prefs=ea_prefs('');

if ~exist(prefs.lc.datadir,'file');
    if ~(isdeployed && ~isfield(prefs,'firstrun')) % dont enter first time in deployed app
        try mkdir(prefs.lc.datadir); end
        if ~exist(prefs.lc.datadir,'file')
            warning('Connectome data directory could not be established. Please set prefs.lc.data dir folder to a valid directory in preferences.');
        end
    end
end
if ~isempty(prefs.ltx.pdfconverter)
   if ~exist(prefs.ltx.pdfconverter,'file')
      warning('LaTeX PDF converter not set properly. Please edit prefs.ltx.pdfconverter in preferences.'); 
   end
end
if ~isempty(prefs.ls.dir)
    if ~exist(prefs.ls.dir,'file')
        ea_mkdir(prefs.ls.dir);
        if ~exist(prefs.ls.dir,'file')
            warning('Lead DBS Server directory could not be established. Please set prefs.ls.dir folder to a valid directory in preferences.');
        end
        
    end
end
if ~isempty(prefs.ixi.dir)
    if ~exist(prefs.ixi.dir,'file')
        warning('IXI data directory empty. Please set prefs.ixi.dir folder to a valid directory in preferences.');
        try
            mkdir(prefs.ixi.dir);
        catch
            warning('IXI data directory could not be established. Please set prefs.ixi.dir folder to a valid directory in preferences.');
        end
    end
end