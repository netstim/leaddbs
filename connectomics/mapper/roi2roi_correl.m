function roi2roi_correl(sublist,roilist,outfile)

% read in text files to create matlab cells
sID=fopen(sublist);
sublist=textscan(sID,'%s');
sublist=sublist{1};
fclose(sID);

rID=fopen(roilist);
roilist=textscan(rID,'%s');
roilist=roilist{1};
fclose(rID);

% iterate subjects/rois to create matrix.
for sub=1:length(sublist)
                [sub_dir,sub_name]=fileparts(sublist{sub});
                
    for roi=1:length(roilist)
            
            tc(:,roi)=load([sub_dir,filesep,sub_name,filesep,'fcMRI_ANALYSIS',filesep,sub_name,'_',roilist{roi},'.dat']);
            
        
    end    
    R(:,:,sub)=corr(tc);
end
meanR=mean(R,3);

Fz=atanh(R);
meanFz=mean(Fz,3);

for xx=1:length(roilist)
   for yy=1:length(roilist) 
    [~,~,~,stats]=ttest(Fz(xx,yy,:));
    T(xx,yy)=stats.tstat;
   end
end

save(outfile,'R','Fz','T','meanFz','meanR');