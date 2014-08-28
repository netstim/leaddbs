function res=ea_update(varargin)
res=(ea_getvsn('local')==ea_getvsn('web'));

if nargin
    
    if strcmp(varargin{1},'force')
       res=[0,0]; 
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
    if ~res(2) % update content
        mkdir([earoot,'tmp'])
        disp('Downloading content, this may take a while...');
        urlwrite('http://www.lead-dbs.org/release/lead_content.zip',[earoot,'tmp',filesep,'lead_content.zip']);
        disp('Extracting content...');
        try
            unzip([earoot,'tmp',filesep,'lead_content.zip'],[earoot,'tmp',filesep]);
            
            delete([earoot,'tmp',filesep,'lead_content.zip']);
            disp('Moving content to the right place...');
            fis=dir([earoot,'tmp',filesep]);
            
            % move template folder..
            movefile([earoot,'tmp',filesep,'templates'],[earoot,fis(fi).name,'templates']);
            % move content from atlases..
            fis=dir([earoot,'tmp',filesep,'atlases',filesep]);
            for fi=1:length(fis)
                if ~strcmp(fis(fi).name,'.') && ~strcmp(fis(fi).name,'..')
                    movefile([earoot,'tmp',filesep,'atlases',filesep,fis(fi).name],[earoot,filesep,'atlases',filesep,fis(fi).name]);
                end
            end
            
            % move content from fibers..
            fis=dir([earoot,'tmp',filesep,'fibers',filesep]);
            for fi=1:length(fis)
                if ~strcmp(fis(fi).name,'.') && ~strcmp(fis(fi).name,'..')
                    movefile([earoot,'tmp',filesep,'fibers',filesep,fis(fi).name],[earoot,filesep,'fibers',filesep,fis(fi).name]);
                end
            end
        catch
            cd(earoot)
            system(['unzip ',earoot,'tmp',filesep,'lead_content.zip']);
        end
        
        disp('Cleaning up...');
        rmdir([earoot,'tmp'],'s')
        disp('Done.');
    end
else
    disp('LEAD-DBS is already up-to-date.');
end




