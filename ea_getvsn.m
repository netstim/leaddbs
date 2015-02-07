function version=ea_getvsn(com)
% This function simply exports the version of the current Lead
% distribution. For updates please see lead-dbs.org
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

switch com
    
    case 'web'
        try
        urlwrite('http://www.lead-dbs.org/release/.version.txt','.webversion.txt','Timeout',5);
        load('.webversion.txt');
        delete('.webversion.txt');
        version=X_webversion;
        catch
            version=[nan;nan];
        end
    case 'local'
        try
            load([fileparts(which('lead')),filesep,'.version.txt']);
            
            version=X_version;
            
        catch % initialize versioning (if .version.txt has gotten lost).
            version=[1.0;1.0];
            fileID = fopen([fileparts(which('lead')),filesep,'.version.txt'],'w');
            fprintf(fileID,'%6.3f\n',version);
            fclose(fileID);
        end
end
version=version(1);


% rough version history: (see git repo for details)
% 
% vstr='v 1.2';
% added DICOM import
% vstr='v 1.15';
% improved and generalized electrode tip recognition, behaviour and added auto-learning feature.
% vstr='v 1.05';
% added lead_group.
% vstr='v 0.91';
% refined manual correction possibilities.
% vstr='v 0.9';
% added support for group-studies (different colors for each group in
% electrode-renderings, fibers and connected anatomical regions).
% vstr='v 0.85';
% added isovalue rendering possibiliy (the hullmethod can be set in
% ea_prefs).
% vstr='v 0.8';
% changed name to Lead DBS (Localization of Electrodes After DBS-Surgery).
% vstr='v 0.7';
% added possibility to display more than one pair of electrodes.
% vstr='v 0.65';
% added Cg25 support
% vstr='v 0.6';
% added implementation for monopolar VAT-reconstruction following the
% approach of maedler 2012
% vstr='v 0.56';
% added normalization routines after Witt 2013 and Schoenecker 2009 that
% use preoperative images.
% vstr='v 0.55';
% added more accurate wireframe reconstruction algorithm for render view as
% implemented in ea_concavehull. possible to use this by changing ea_prefs
% vstr='v 0.5';
% added CT support
% vstr='v 0.48';
% added normalization routine after Schoenecker 2009
% vstr='v 0.47';
% added GUI and possibility to save each step of
% reconstruction/visualization
% vstr='v 0.46';
% added possibility to visualize fiber connectivities
% vstr='v 0.45';
% added possibility to visualize fibers
% vstr='v 0.4';
% added ability of render view
% vstr='v 0.35';
% added electrode tip templates to largely improve ability for height
% reconstruction
% vstr='v 0.33';
% added possibility to change volume contrasts.
% vstr='v 0.3';
% added ability to export slices, large robustness improvements
% vstr='v 0.25';
% first manual refinement possible
% vstr='v 0.2';
% rewrote whole program based on first experiences
% vstr='v 0.1';
% very raw version of automatic reconstruction of electrodes trajectories
