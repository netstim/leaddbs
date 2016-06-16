  
      
function data = loadmask(fn,threshold,ref,edges)
        data = [];
        
        if iscell(fn),
            [fp fndum fext] = fileparts(fn{1});
        else
            [fp fndum fext] = fileparts(fn);
        end;
        if strcmp(fext(1:4),'.nii') || strcmp(fext(1:4),'.hdr'),       
           if size(edges,1) == 4, % a edges matrix, 
               mrdum.dataAy = ref;
               mrdum.edges = edges;
           else % a nifti header
               mrdum =edges.hdr;
           end;
                     
           if iscell(fn)
               masknii = load_untouch_nii(fn{2});
               maskdata = reslice_nifti(masknii.img,masknii.hdr,mrdum,2);
               masknii = load_untouch_nii(fn{1});
               maskdata(:,:,:,2) = reslice_nifti(masknii.img,masknii.hdr,mrdum,2);
           else
               masknii = load_untouch_nii(fn);
               maskdata = reslice_nifti(masknii.img,masknii.hdr,mrdum,2);
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
            
        else
            if not(iscell(fn)),
                fn = {fn};
            else
                if length(fn)>1,
                    fn = fn([2 1]);
                end;
            end;
            
            for k = 1:length(fn),
                mrstruct = load(fn{k});
                if isfield(mrstruct,'maskNamesCell')
                    mastr = mrstruct;
                    if isempty(threshold),
                        [res ok] = listdlg('ListString',mastr.maskNamesCell,'SelectionMode','single','Name','Select a Mask');
                    else
                        ok = true;
                        res = threshold;
                    end;
                    if ok,                    
                        maskdata(:,:,:,k) = mastr.maskCell{res};                    
                    end;
                else                            
                    if not(isfield(mrstruct,'mrStruct')),
                        errordlg('Error while reading','Loading Mask');
                        return;
                    else
                       mrstruct = mrstruct.mrStruct;
                       if length(threshold) < k,
                           if not(islogical(mrstruct.dataAy)),
                                threshold_tmp = chooseThreshold_stackview(ref,mrstruct.dataAy(:,:,:,1));
                                if isempty(threshold_tmp)
                                   return;                          
                                end;                           
                            else
                                threshold_tmp = 0.5;
                            end; 
                            threshold(k) = threshold_tmp;
                       end;
                       maskdata(:,:,:,k) = mrstruct.dataAy>threshold(k);
                    end;
                end;
            end;
        end;
        
        data.mask = maskdata;
        data.threshold = threshold;
        
        
       