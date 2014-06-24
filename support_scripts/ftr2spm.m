% put FTR-File to mm and spm notation.
clc
clear
nii=load_nii('/mnt/marc/HornAndreas/MATLAB/fspm8/templates/T1.nii');
ysize=size(nii.img,2)+1;
clear nii

ftr=load('/mnt/marc/HornAndreas/topography/ggomprFTR_20fold.mat');



gibbsconnectome=ftr.curveSegCell;

for fib=1:length(ftr.curveSegCell)
    
    gibbsconnectome{fib}=[ftr.curveSegCell{fib}(:,2),ysize-ftr.curveSegCell{fib}(:,1),ftr.curveSegCell{fib}(:,3),ones(size(ftr.curveSegCell{fib},1),1)];

    gibbsconnectome{fib}=map_coords(gibbsconnectome{fib}','/mnt/marc/HornAndreas/MATLAB/fspm8/templates/T1.nii')';
    
    
end

save('gibbsconnectome.mat','gibbsconnectome');