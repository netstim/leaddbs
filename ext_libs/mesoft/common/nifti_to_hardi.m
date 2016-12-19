%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%   hardi = nifti_to_hardi(input,T)
%
%   converts a series of niftis into a mrstruct (HARDI)
%
%   input - 
%      a) a txt filename containing filenames and 
%      gradient information as follows:
%          file1.nii dir1x dir1y dir1z b1
%          file2.nii dir2x dir2y dir2z b2
%            .....
%
%      b) a cellarray of strings containing
%      nifti-filenames, where an additionally .mat is expected to exist
%      which contains the gradient info. It should has to have the same 
%      name like the nifti and contain a structure with a field userdata.b 
%      and userdata.g
%
%    T - optional, containing a 3x3 matrix which is applied onto the 
%     gradient directions,
%
%    hardi - output containing the hardi structure, e.g. save it with 
%    mrstruct_write(hardi,'mynew_HARDI.mat');
%


function hardi = nifti_to_hardi(input,T)

if iscell(input),
   hardi = nifti_to_hardi_from_immat(input);
elseif isstr(input),
   hardi =  nifti_to_hardi_from_imcom(input);
else
    error('nifti_to_hardi: wrong input format');
    return;
end;

 
if not(exist('T'))
    T = 1;
end;



for k = 1:size(hardi.user.bDir,2),
    hardi.user.bDir(:,k) = T*hardi.user.bDir(:,k);
    hardi.user.bTensor(:,:,k) = T*hardi.user.bTensor(:,:,k)*T';
end;
    


function hardi = nifti_to_hardi_from_immat(niftifis)

for k = 1:length(niftifis),
    [p n e] = fileparts(niftifis{k});
    info = load(fullfile(p,[n '.mat']));
    
    [ U,~,V] = svd(info.userdata.mat(1:3,1:3));
    M = [0 -1 0; 1 0 0 ; 0 0 1]*(V*U')*info.userdata.ref(1:3,1:3);
    
    bval(k) = info.userdata.b;
    dir(:,k) = M*info.userdata.g';
end;




hardi = n_to_h(niftifis,dir,bval);

    

function hardi = nifti_to_hardi_from_imcom(imagecom_fi)

[p n e] = fileparts(imagecom_fi);
fi = fopen(imagecom_fi,'r');
c = textscan(fi,'%s %s %s %s %s');

for k = 1:length(c{1});
    names{k} = fullfile(p,c{1}{k});
    bval(k) = sscanf(c{2}{k},'b=%i');
    if bval(k) > 0,
        d1 = str2num(c{3}{k});
        d2 = str2num(c{4}{k});
        d3 = str2num(c{5}{k});
        dir(:,k) = [d1 d2 d3]';
    else
        dir(:,k) = [0 0 0]';
    end;
end;

hardi = n_to_h(names,dir,bval);
    

return;


function hardi = n_to_h(data_names,bdirs,bvals);

[hardi err] = nifti_to_mrstruct('series3D',data_names);
if not(isempty(err))
    warning(err);
end;

hardi.user.bfactor = bvals(:);

%bdirs = bdirs([2 1 3],:);
bdirs = bdirs./repmat(eps+sqrt(sum(bdirs.^2)),[3 1]);

hardi.user.bDir = bdirs;
for k = 1:size(hardi.user.bDir,2),
    hardi.user.bTensor(:,:,k) = hardi.user.bDir(:,k)*hardi.user.bDir(:,k)' *bvals(k);
end;

return;

