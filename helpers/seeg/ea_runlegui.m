function run = ea_runlegui(options)
	prefsPath = options.bids.getPrefs(options.bids.subjId{1}, 'uiprefs', 'mat');
    prefs = load(prefsPath);
    prefs.reconmethod = 'LeGUI (Davis 2021)';
    save(prefsPath, '-struct', 'prefs');
    open_legui(fullfile(options.root, options.patientname));
end

function open_legui(selDir)
    arguments
        selDir (1,1) string
    end
    if ~isfolder(selDir), error("Not a folder: %s", selDir); end

    % 1) Start LeGUI
    app = LeGUI(selDir);
    %     tmpDir = tempname; mkdir(tmpDir);
%     shadowFile = fullfile(tmpDir, 'LeG_lastDir.m');
%     fid = fopen(shadowFile, 'w'); assert(fid>0, 'Could not create shadow LeG_lastDir.m');
%     fprintf(fid, 'function p = LeG_lastDir(varargin)\n');
%     fprintf(fid, '%% Shadowed picker: always return provided path\n');
%     fprintf(fid, 'p = ''%s'';\n', strrep(selDir, '''', '''''')); % escape quotes
%     fprintf(fid, 'end\n');
%     fclose(fid);
% 
%     addpath(tmpDir, '-begin');
%     c = onCleanup(@()cleanupShadow(tmpDir));
%     cb = app.LoadImgsBtnH.ButtonPushedFcn;
%     if isempty(cb), error('Could not find Load Images button callback on the app.'); end
%     if iscell(cb) && isa(cb{1}, 'function_handle')
%         feval(cb{1}, app, []); 
%     else
%         feval(cb, app, []);
%     end
end


function cleanupShadow(tmpDir)
    % Remove temp path and delete the temp folder
    if ~isempty(tmpDir) && isfolder(tmpDir)
        rmpath(tmpDir);
        try
            delete(fullfile(tmpDir, 'LeG_lastDir.m'));
            rmdir(tmpDir);
        catch
            % best-effort cleanup
        end
    end
end
