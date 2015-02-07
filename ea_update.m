function res=ea_update(varargin)
res=(ea_getvsn('local')==ea_getvsn('web'));



if nargin
    
    if strcmp(varargin{1},'force')
       res=0; 
    else
    return
    end
end

earoot=[fileparts(which('lead')),filesep];

if any(~res)
    if ~res(1) % update code
        mkdir([earoot,'tmp'])
        disp('Downloading code...');
        urlwrite('http://www.lead-dbs.org/release/lead_dbs.zip',[earoot,'tmp',filesep,'lead_dbs.zip']);
        disp('Extracting code...');
        try
            unzip([earoot,'tmp',filesep,'lead_dbs.zip'],[earoot,'tmp',filesep]);
            delete([earoot,'tmp',filesep,'lead_dbs.zip']);
            disp('Moving code to the right place...');
            fis=dir([earoot,'tmp',filesep]);
            for fi=1:length(fis)
                if ~fis(fi).isdir || strcmp(fis(fi).name,'ext_libs') || strcmp(fis(fi).name,'icons') || strcmp(fis(fi).name,'nii_lib') || strcmp(fis(fi).name,'support_scripts')
                    movefile([earoot,'tmp',filesep,fis(fi).name],[earoot,fis(fi).name]);
                end
            end
        catch
            cd(earoot)
            system(['unzip ',earoot,'tmp',filesep,'lead_dbs.zip']);
        end
        disp('Cleaning up...');
        rmdir([earoot,'tmp'],'s')
        disp('Done.');
    end
    disp('Writing version file.');
    version=ea_getvsn('web');
    if ~isnan(version)
    fileID = fopen([fileparts(which('lead')),filesep,'.version.txt'],'w');
    fprintf(fileID,'%6.3f\n',version);
    fclose(fileID);
    end
    disp('Restarting LEAD-DBS.');
    lead;
else
    disp('LEAD-DBS is already up-to-date.');
end




