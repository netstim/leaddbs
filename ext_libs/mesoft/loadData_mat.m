function data = loadData_mat(dwifile,gradfile_dummy,maskfile,threshold)

       data = [];
       if not(exist('dwifile')),
            [files path] = uigetfile({'*_HARDI.mat', 'All accepted filetypes';'*.*','All files'},'Select HARDI data to open');      
            if files==0,
                return;
            end;
            dwifile = fullfile(path,files);
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
                
       mrstruct = load(dwifile);
       if not(isfield(mrstruct,'mrStruct')),
          errordlg('Error while reading','Loading DWI');                        
          return
       end
       mrstruct = mrstruct.mrStruct;
          
       if ~isfield(mrstruct.user,'bTensor'),
           ready;
           errordlg('Gradient direction data missing!','Loading MAT data');
           return;
       end;
       
       if max(mrstruct.user.bTensor(1,1,:)+mrstruct.user.bTensor(2,2,:)+mrstruct.user.bTensor(3,3,:)) < 1.001,
           mrstruct.user.bTensor = mrstruct.user.bTensor *mrstruct.user.bfactor;
       end;
       
       

       data.dwifile = dwifile;
       data.dwi = mrstruct.dataAy;
       data.gradfile = [];
       data.tensor = mrstruct.user.bTensor;
       data.edges = mrstruct.edges;
       data.vox = mrstruct.vox;
       data.name = mrstruct.patient;
       
       
       if not(exist('threshold')),
           threshold = [];
       end;
       
       data.WM = loadmask(maskfile,threshold,data.dwi(:,:,:,1),mrstruct.edges);
       data.WM.file = maskfile;
       
       
     