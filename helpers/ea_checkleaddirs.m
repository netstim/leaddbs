function ea_checkleaddirs
prefs=ea_prefs('');
try
if ~exist(prefs.lc.datadir,'file');
    mkdir(prefs.lc.datadir);
    if ~exist(prefs.lc.datadir,'file')
        ea_warning('Connectome data directory could not be established. Please set prefs.lc.datadir folder to a valid directory in preferences.');
    end
end
end
try
   if ~exist(prefs.ltx.pdfconverter,'file')
      warning('LaTeX PDF converter not set properly. Please edit prefs.ltx.pdfconverter in preferences.'); 
   end
end
try
    if ~exist(prefs.ls.dir,'file')
        mkdir(prefs.ls.dir);
        if ~exist(prefs.ls.dir,'file')
            ea_warning('Lead DBS Server directory could not be established. Please set prefs.ls.dir folder to a valid directory in preferences.');
        end
        
    end
end
try
    if ~exist(prefs.ixi.dir,'file')
        warning('IXI data directory empty. Please set prefs.ixi.dir folder to a valid directory in preferences.');
        mkdir(prefs.ixi.dir);
        if ~exist(prefs.ixi.dir,'file')
            warning('IXI data directory could not be established. Please set prefs.ixi.dir folder to a valid directory in preferences.');
        end
    end
end