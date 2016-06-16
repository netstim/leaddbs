function data = loadData_nii(dwifile,gradfile,maskfile,threshold)
    
       data = [];

       if not(exist('dwifile')),
            [files path] = uigetfile({'*.nii;*.nii.gz;*.hdr', 'All accepted filetypes'},'Select nii containing DWI information','MultiSelect', 'on');      
            if not(iscell(files)),
                if files==0,
                    return;
                end;
            end;
            if iscell(files),
                dwifile = cellfun(@(x) fullfile(path,x),files,'uniformoutput',false);
            else              
                dwifile = fullfile(path,files);
            end;
            
            [gradf pathgradf] = uigetfile({'*', 'All accepted filetypes'},'Select gradient information files','MultiSelect', 'on');      
            if length(gradf)~=2,
                errordlg('You have to select 2 files containing gradient directions and b-values (FSL style)');
                return;
            end;
            gradfile{1} = fullfile(pathgradf,gradf{1});
            gradfile{2} = fullfile(pathgradf,gradf{2});
            
            
            [fnmask pathmask] = uigetfile({'*.mat;*.nii','Accepted Files (*.mat,*.nii,*.hdr)'},'Load mrStruct/maskStruct','MultiSelect','on');
            if isempty(fnmask),
                 return;
            end;
            if iscell(fnmask),
                maskfile{1} = fullfile(pathmask,fnmask{1});
                maskfile{2} = fullfile(pathmask,fnmask{2});
            else
                maskfile = fullfile(pathmask,fnmask);
            end;
            
       end;
                
       if iscell(dwifile),            
            dwinii = load_untouch_nii(dwifile{1});
            dwinii.img = zeros([size(dwinii.img) length(dwifile)]);
            for k = 1:length(dwifile),
                tmp = load_untouch_nii(dwifile{k});
                dwinii.img(:,:,:,k) = tmp.img;
            end;                
       else
            dwinii = load_untouch_nii(dwifile);
       end;
       
       gradf1 = importdata(gradfile{1});
       gradf2 = importdata(gradfile{2});
       if size(gradf1,1) == 3,
           bvec = gradf1;
           bval = gradf2;
       else
           bvec = gradf2;
           bval = gradf1;
       end;

       T = 1;
              
       for k = 1:size(bval,2),
           gdir = T*bvec(:,k);
           gdir = gdir / (eps+norm(gdir));
           tensor(:,:,k) = gdir*gdir' *bval(k);
       end;
       
       
       data.dwifile = dwifile;
       data.dwi = dwinii.img;
       data.gradfile = gradfile;
       data.tensor = tensor;
       data.edges = rmfield(dwinii,'img');      
       data.vox = dwinii.hdr.dime.pixdim(2:4);
       data.name = dwinii.fileprefix;
       
              
   
       if not(exist('threshold')),
           threshold = [];
       end;
              
       
       data.WM = loadmask(maskfile,threshold,data.dwi(:,:,:,1),data.edges.hdr);
       data.WM.file = maskfile;
       
       
       
      
function data = loadmask(fn,threshold,ref,edges)
        data = [];
        
        if iscell(fn),
            [fp fndum fext] = fileparts(fn{1});
        else
            [fp fndum fext] = fileparts(fn);
        end;
        if strcmp(fext(1:min(4,end)),'.gz') ||  strcmp(fext(1:min(4,end)),'.nii') || strcmp(fext(1:min(4,end)),'.hdr'),                     
                     
           if iscell(fn)
               masknii = load_untouch_nii(fn{2});
               maskdata = reslice_nifti(masknii.img,masknii.hdr,edges,2);
               masknii = load_untouch_nii(fn{1});
               maskdata(:,:,:,2) = reslice_nifti(masknii.img,masknii.hdr,edges,2);
           else
               masknii = load_untouch_nii(fn);
               maskdata = reslice_nifti(masknii.img,masknii.hdr,edges,2);
           end;
           for k = 1:size(maskdata,4),
               if length(threshold) < k,
                    threshold_tmp = chooseThreshold_stackview(ref,maskdata(:,:,:,k));
                    if isempty(threshold_tmp)
                       return;
                    end;                    
                    threshold(k) = threshold_tmp;
               end;
               maskdata(:,:,:,k) = maskdata(:,:,:,k) > threshold(k);
           end;                
            
        end;
        
        data.mask = maskdata;
        data.threshold = threshold;
        
        
