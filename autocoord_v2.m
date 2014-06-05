function results=autocoord_v2(patientname,options)


% This function is the main function of eAuto-DBS. It will generate a 8x3
% vector of coordinates. Rows 1-4 will be the right electrodes, rows 5-8
% the left ones. fitline{1} will be the right trajectory, fitline{2} the
% left one.



%%

if options.normalize
   eauto_normalize(options.root,patientname); 
    
end

if ~exist([options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii'],'file')
    try
        copyfile([options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii.gz'],[options.root,patientname,filesep,patientname,'_tra_brain_A3_final2.nii.gz']);
    end
    system(['gunzip ',options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii.gz']);
    try
        movefile([options.root,patientname,filesep,patientname,'_tra_brain_A3_final2.nii.gz'],[options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii.gz']);
    end
    
end
if ~exist([options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt.nii'],'file')
    try
        copyfile([options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt.nii.gz'],[options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt2.nii.gz']);
    end
    system(['gunzip ',options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt.nii.gz']);
    try
        movefile([options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt2.nii.gz'],[options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt.nii.gz']);
    end
    
end

for side=options.sides
%try    
    % call main routine reconstructing trajectory for one side.
    [coords,trajvector{side},fitline{side}]=autocoord_side(patientname,options,side);
    
    
    
    %% refit electrodes starting from first electrode (this is redundant at
    %% this point.
    
    
    
    coords_mm = map_coords(coords', [options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii'])';
    % rename matfile to text
    movefile([options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.mat'],[options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.txt']);
    % load matrices
    tramat=load([options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.txt']);
    % rename the file to .mat again
    movefile([options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.txt'],[options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.mat']);
    
    
    [d,distmm]=calc_distance(options.eldist,trajvector{side},tramat(1:3,1:3),[options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii']);
    
    comp = map_coords([0,0,0;trajvector{side}]', [options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii'])'; % (XYZ_mm unaltered)
    
    trajvector{side}=diff(comp);
    
    
    
    normtrajvector{side}=trajvector{side}./norm(trajvector{side});
    
    for electrode=2:4
        
        coords_mm(electrode,:)=coords_mm(1,:)-normtrajvector{side}.*((electrode-1)*distmm);
        
    end
    
    
    
   
    
%end
 if side==1

     if ~exist('coords_mm','var')
         coords_mm=nan(4,3);
     end
     rightcoords_mm=coords_mm;
    elseif side==2
        if ~exist('coords_mm','var')
            coords_mm=nan(4,3);
        end
        
            coords_mm=[rightcoords_mm;coords_mm];       % append right electrode contacts to left ones.
        
    end
end



if length(coords_mm)==4 % only one side was processed.
    if options.sides==1
        coords_mm=[nan(4,3);coords_mm];
    elseif options.sides==2
        coords_mm=[coords_mm;nan(4,3)];
    end
end





try
    %realcoords=load([options.root,patientname,filesep,'L.csv']);
    %realcoords=realcoords(:,1:3);
    
    realcoords=read_fiducials([options.root,patientname,filesep]);
catch
    showdis('No manual survey available.',options.verbose);
realcoords=[];
end


[resultfig,coords_mm]=showresultfig(coords_mm,realcoords,fitline,[options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt.nii'],patientname,options);



%% check traject sanity

for side=options.sides
try
    trajectissane=checktrajectsanity(trajvector{side});
    if ~trajectissane
       disp(['Trajectory of side ',num2str(side),' seems not to have been correctly reconstructed. Check manually.']); 
    end

end
end




if options.verbose>1; saveas(resultfig,[options.root,patientname,filesep,'automan.png']); end
if options.verbose>1; saveas(resultfig,[options.root,patientname,filesep,'automan.fig']); end

if options.verbose>3; close(resultfig); end

results.coords_mm=coords_mm;
try
    results.realcoords=realcoords;
    
    
    
    for electrode=1:length(coords_mm)
        
        results.distances(electrode)=pdist([coords_mm(electrode,:);realcoords(electrode,:)]);
        
    end
    results.fit=nanmean(results.distances);
    
end









