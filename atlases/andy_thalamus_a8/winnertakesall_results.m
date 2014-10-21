%% analysis to map thalamus connectivities..
function group_results
clc
root='/PA/Neuro/_projects/lead/lead_dbs/atlases/andy_thalamus_a8/mixed/';
spmroot=[spm('dir'),filesep];

export_normalized=0;
export_zscored=1;
comps={'frontotemporal','motor','parietooccipital','sensory'};
%comps={'occipital','postparietal','prefrontal','premotor','primarymotor','sensory'};

%comps={'postparietal','prefrontal','premotor','primarymotor','sensory','temporal'};

for fi=1:length(comps)
   fis{fi}=[root,'scmean_',comps{fi},'.nii,1'];
end

matlabbatch{1}.spm.util.cat.vols = fis;
matlabbatch{1}.spm.util.cat.name = 'compound.nii';
matlabbatch{1}.spm.util.cat.dtype = 4;

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear jobs matlabbatch

  V=spm_vol([root,'compound.nii']);
   X=spm_read_vols(V);
   Y=zeros(size(X));
for xx=1:size(X,1)
    for yy=1:size(X,2)
        for zz=1:size(X,3)
            if any(X(xx,yy,zz,:))
            [~,whowins]=max(X(xx,yy,zz,:));
            Y(xx,yy,zz,whowins)=1;
            end
        end
    end
end

for vol=1:length(V);
   V(vol).fname=[root,'compound_wta.nii']; 
end
nii=load_untouch_nii([root,'compound.nii']);
nii.img=Y;
save_untouch_nii(nii,[root,'compound_wta.nii']);
spm_file_split([root,'compound_wta.nii']);


for fi=1:length(comps)
movefile([root,'compound_wta_0000',num2str(fi),'.nii'],[root,'compound_thalamus_',comps{fi},'.nii']);
end


function  th_dispercent(varargin)
%
percent=round(varargin{1}*100);

if nargin==2
    if strcmp(varargin{2},'end')
        fprintf('\n')
        fprintf('\n')
        
        fprintf('\n')
        
    else
        fprintf(1,[varargin{2},':     ']);
        
        
    end
else
    fprintf(1,[repmat('\b',1,(length(num2str(percent))+1)),'%d','%%'],percent);
end


