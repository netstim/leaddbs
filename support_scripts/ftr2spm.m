% put FTR-File (from Global Tracker of DTI&Fibertools for SPM) to mm and spm notation.
clc
clear
% load corresponding anatomical image (could be the b0 image)
b0file=[fileparts(which('spm')),filesep,'templates',filesep,'T1.nii'];
nii=load_nii(b0file);
ysize=size(nii.img,2)+1;
clear nii

thinnout=[5,10,50,100,500]; % leave empty if you don't need lightweight copies.

% load ftr-file
ftr=load('ggomprFTR_20fold.mat');



gibbsconnectome=ftr.curveSegCell;

for fib=1:length(ftr.curveSegCell)
    
    gibbsconnectome{fib}=[ftr.curveSegCell{fib}(:,2),ysize-ftr.curveSegCell{fib}(:,1),ftr.curveSegCell{fib}(:,3),ones(size(ftr.curveSegCell{fib},1),1)];

    gibbsconnectome{fib}=map_coords(gibbsconnectome{fib}',b0file)';
    gibbsconnectome{fib}=gibbsconnectome{fib}(:,1:3);
    
    
end

% save fiberset.
save('gibbsconnectome.mat','gibbsconnectome','-v7.3');


% if you want to export light-weight copies of the fiberset, set thinning
% to a vector determining how many fibers to skip for each fiber to leave
% in.
if ~isempty(thinnout)
    gc=gibbsconnectome;
    for t=1:length(thinnout)
        
   gibbsconnectome=gc{1:thinnout(t):end};
   save(['gibbsconnectome_',num2str(thinnout(t)),'.mat'],'gibbsconnectome');
        
        
        
    end
    
end