function ea_robustaverage_nii(filelist,outputname,components)

if ~exist('components','var')
    components=3;
    % by default assumes niis with 3 components (i.e. X x Y x Z x 3 dim files).
end

% load nii headers if necessary
for d=1:length(filelist)
    for i=1:components
    V{d}(i)=ea_open_vol([filelist{d},',1,',num2str(i)]);
    end
end
O = V{1};
for i=1:length(O)
    O(i).fname=[outputname];
end

O = spm_create_vol(O);
ea_dispercent(0,'Robust averaging nifti files');
for i=1:numel(V{1});
    for p=1:V{1}(1).dim(3)
        for d=1:length(V);
            Y(:,:,d) = spm_slice_vol(V{d}(i),spm_matrix([0 0 p]),V{d}(i).dim(1:2),0);
        end
        Y=ea_robustmean(Y,3);
        %Y=mean(Y,3);
        O(i) = spm_write_plane(O(i),Y,p);
        ea_dispercent(((p/V{1}(1).dim(3)))/numel(V{1}));
    end
end
ea_dispercent(1,'end');