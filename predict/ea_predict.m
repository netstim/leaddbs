function ea_predict(options)
% main function called from lead predict

funct=options.predict.model_mfile;
ea_dispercent(0,'Predicting patient outcome');
for pt=1:length(options.uivatdirs)
   
    feval(funct,options,pt);
   ea_dispercent(pt/length(options.uivatdirs)); 
end
ea_dispercent(1,'end');



