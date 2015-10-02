function lutEntry = createInteractionLUTs(ten,sinterp,Pstruc,fixedtissuemodel,thedir)

verbose = false;



% model model correlations
fid = fopen(fullfile(thedir,'interactionLUTs.h'),'w');
fprintf(fid,'/********* DWscheme dependent model-model interaction ***********/\n\n');
bestpp = createModelModelLUT(ten,Pstruc,fixedtissuemodel,fid,verbose);

% model signal correlations
fprintf(fid,'\n\n/********* DWscheme dependent model-signal interaction ***********/\n\n');
model2alphaLUT = createModelSignalLUT(ten,Pstruc,fixedtissuemodel,sinterp,fid,verbose);


fclose(fid);


lutEntry.model2alphaLUT = model2alphaLUT;




    
    
    