function ea_dicom_import(options)
% This function converts DICOM files in your in directory and outputs them
% to the out directory as specified by lead. This function is pretty much
% specialized to the DICOM format as used @ Charite University Medicine and
% might not work as expected in your clinical setting. You can
% alternatively import DICOM files using software like SPM or DCM2NII.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

disp('Importing DICOM files...');

global dcfilename
global tmpoutdir
indir=options.prefs.lp.dicom.infolder;
outdir=options.prefs.lp.dicom.outfolder;

f=dir(indir);
% zipfile support..
for scan=1:length(f)
   [p,fi,e]=fileparts(f(scan).name);
   if strcmp(e,'.zip')
       unzip([indir,f(scan).name],indir);
       delete([indir,f(scan).name]);
   end
end



f=dir(indir);




for scan=1:length(f)
    if f(scan).isdir && ~strcmp(f(scan).name,'.') && ~strcmp(f(scan).name,'..')
        
        namecell=textscan(f(scan).name,'%s','Delimiter',',');
        
        nachname=char(namecell{1}(1));
        try
            vorname=char(namecell{1}(2));
            
            
            vornamecell=textscan(vorname,'%s','Delimiter','_');
            
            vorname=char(vornamecell{1}(1));
        catch
            vorname='';
        end
        
        
        name=[nachname,vorname];
        
        lpath=genpath([indir,f(scan).name]);
        delims=strfind(lpath,':');
        from=1;
        for d=1:length(delims)
            pfolds{d}=lpath(from:delims(d)-1);
            from=delims(d)+1;
        end
        
        
        %% dicom import
        
        % create output folder in working directory..
        if exist([outdir,name],'file')==false
            mkdir([outdir,name]);
        end
        tmpoutdir=[outdir,name];
        
        fcnt=1;
        for pfold=1:length(pfolds)
            files=dir(pfolds{pfold});
            for i = 1:(length(files))
                if ~strcmp(files(i).name(1),'.') && ~files(i).isdir
                    filecell{fcnt}=[pfolds{pfold},filesep,files(i).name];
                    fcnt=fcnt+1;
                end
            end
        end
      
        
        matlabbatch{1}.spm.util.dicom.data = filecell;
        %%
        clear filecell
        
        matlabbatch{1}.spm.util.dicom.root = 'flat';
        matlabbatch{1}.spm.util.dicom.outdir = {tmpoutdir};
        matlabbatch{1}.spm.util.dicom.convopts.format = 'nii';
        matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;
        
        jobs(1)=matlabbatch;
          
        try
            cfg_util('run',jobs);
            worked=1;
        catch
            warning('Invalid DICOM structure.');
            worked=0;
        end
        clear matlabbatch jobs
        
        
        if worked
            % cleanup DICOMs first
            if options.prefs.dicom.dicomfiles==1 % delete DICOMs
                rmdir([indir,f(scan).name],'s');
            elseif options.prefs.dicom.dicomfiles==2 % move DICOMs
                mkdir([tmpoutdir,filesep,'DICOM']);
                movefile([indir,f(scan).name],[tmpoutdir,filesep,'DICOM',filesep,f(scan).name]);
            end
            
            
            % display file:
            nid=dir([tmpoutdir,filesep,'s*.nii']);
            for ni=1:length(nid)
                if ~strcmp(nid(ni).name,'sag.nii')
                    movefile([tmpoutdir,filesep,nid(ni).name],[tmpoutdir,filesep,'t',nid(ni).name]) % remove from search key (s*)..
                    dcfilename=[tmpoutdir,filesep,'t',nid(ni).name];
                    ea_imageclassifier;
                    
                end
            end
            
            
            
            
        end
    end
    
end