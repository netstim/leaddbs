function ea_dicom_import(options)

global dcfilename
global tmpoutdir
indir=options.prefs.lp.dicom.infolder;
outdir=options.prefs.lp.dicom.outfolder;

f=dir(indir);

for scan=4:length(f)
    
    if f(scan).isdir
        
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
        
        
        %% dicom import
        
        if exist([outdir,name],'file')==false
            mkdir([outdir,name]);
        end
        tmpoutdir=[outdir,name];
        
        
        files=dir([indir,f(scan).name]);
        filecell=cell(length(files)-4,1);
        for i = 4:(length(files)-1)
            filecell{i-3}=files(i).name;
            
        end
        
        single_s_files=cellfun(@(x) [indir,f(scan).name,filesep,x],filecell(:),'Uniformoutput',false);
        
        
        matlabbatch{1}.spm.util.dicom.data = single_s_files;
        %%
        matlabbatch{1}.spm.util.dicom.root = 'flat';
        matlabbatch{1}.spm.util.dicom.outdir = {tmpoutdir};
        matlabbatch{1}.spm.util.dicom.convopts.format = 'nii';
        matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;
        
        jobs(1)=matlabbatch;
        try
            cfg_util('run',jobs);
            worked=1;
        catch
            warning('Invalid DICOM folder.');
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