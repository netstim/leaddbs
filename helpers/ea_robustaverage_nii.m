function ea_robustaverage_nii(filelist,outputname,components)

if ~exist('components','var')
    components=3;
    % by default assumes niis with 3 components (i.e. X x Y x Z x 3 dim files).
end

% load nii headers if necessary
for d=1:length(filelist)
    for i=1:components
    V{d}(i)=spm_vol([filelist{d},',1,',num2str(i)]);
    end
end

usefit=0;
[pth,fn,ext]=fileparts(outputname);
for i=1:3
    O{i} = V{1}(1);
    O{i}.fname=fullfile(pth,[fn,'_',num2str(i),'.nii']);
end



for i=1:3;
    ea_dispercent(0,['Robust averaging dimension ',num2str(i)]);
    O{i} = spm_create_vol(O{i});
    for p=1:V{1}(1).dim(3)
        for d=1:length(V);
            Y(:,:,d) = spm_slice_vol(V{d}(i),spm_matrix([0 0 p]),V{d}(i).dim(1:2),0);
        end
        
        if usefit
            Z=zeros(size(Y,1),size(Y,2));
            
            for xx=1:size(Y,1)
                for yy=1:size(Y,2)
                    Z(xx,yy)=ea_robustmean(squeeze(Y(xx,yy,:)),[],1,1); % potentially decrease k (=5) to make more robust.
                end
            end
        Y=Z;
        else
            Y=ea_robustmean(Y,3,1,0);
        end
        %Y=mean(Y,3);
        O{i} = spm_write_plane(O{i},Y,p);
        ea_dispercent((p/V{1}(1).dim(3)));
    end
    ea_dispercent(1,'end');
end
