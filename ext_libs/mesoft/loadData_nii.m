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
            
            [fnmask pathmask] = uigetfile({'*.nii;*.nii.gz;*.hdr','All accepted filetypes'},'Select nii containing WM mask');
            if fnmask == 0,
                 return;
            end;
            maskfile = fullfile(pathmask,fnmask);
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
       masknii = load_untouch_nii(maskfile);
       
       gradf1 = importdata(gradfile{1});
       gradf2 = importdata(gradfile{2});
       if size(gradf1,1) == 3,
           bvec = gradf1;
           bval = gradf2;
       else
           bvec = gradf2;
           bval = gradf1;
       end;
               
       T = diag([-1 -1 1])*[dwinii.hdr.hist.srow_x(1:3) ; dwinii.hdr.hist.srow_y(1:3) ; dwinii.hdr.hist.srow_z(1:3)];
       [U S V] = svd(T);
       T = U*V';
              
       for k = 1:size(bval,2),
           gdir = T*bvec(:,k);
           gdir = gdir / (eps+norm(gdir));
           tensor(:,:,k) = gdir*gdir' *bval(k);
       end;
       
       res = reslice_nifti(masknii.img,masknii.hdr,dwinii.hdr,2);
       
       if not(exist('threshold'))           
           threshold = chooseThreshold_stackview(dwinii.img(:,:,:,1),res);
       end;
       if isempty(threshold)
           threshold = chooseThreshold_stackview(dwinii.img(:,:,:,1),res);
       end;
       
       data.dwifile = dwifile;
       data.dwi = dwinii.img;
       data.gradfile = gradfile;
       data.tensor = tensor;
       data.edges = rmfield(dwinii,'img');      
       data.vox = dwinii.hdr.dime.pixdim(2:4);
       data.name = dwinii.fileprefix;
       
              
       data.WM.mask = res;
       data.WM.threshold = threshold;
       data.WM.file = maskfile;
       

