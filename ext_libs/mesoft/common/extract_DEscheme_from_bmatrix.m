[file,path] = uigetfile('_info.mat','Select the *_info.mat file of your data set');
load(fullfile(path,file));
DE_scheme = user.PVM_DwDir;
b_tensor = bmatrix/bValue;
DE_size = size(DE_scheme);
new_scheme = zeros([DE_size(1)+2 DE_size(2)]);
B_tensor =zeros(3,3,length(DE_scheme));
for i = 1:DE_size(1)
    tmp = DE_scheme(i,:)' * DE_scheme(i,:);
    B_tensor(:,:,i+2)= b_tensor(:,:,i+2)*trace(tmp);
end
new_scheme = [sqrt(B_tensor(1,1,:)) sqrt(B_tensor(2,2,:)) sqrt(B_tensor(3,3,:))];
new_scheme = new_scheme/max(new_scheme(:));
new_scheme = squeeze(new_scheme);
new_scheme = new_scheme';
save DE_dirs_recalc_from_matrix.txt -ascii new_scheme
