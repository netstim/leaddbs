function p = editParamStruct(p)





    fignum = figure('NumberTitle','off',...
                   'Menubar','none',...
                   'Toolbar','none',...
                   'Name', 'Edit Parameters'); 
    clf;

    pos = get(fignum,'position');
                    

     fnames = fieldnames(p);

     thei = -[0 25 0 0];
     twid = [200 0 0 0];
     
     set(fignum,'position',[pos(1:2) 400 length(fnames)*25 + 80]);
     posFig=get(fignum,'OuterPosition');
     
     pos = [5 length(fnames)*25+40 180 20];
     
     for k = 1:length(fnames),      
         if not(isa(getfield(p,fnames{k}),'function_handle')),
             uicontrol('Style','text','String',fnames{k},'Position',pos,'Tag','fiberGT_text');
             uicontrol('Style','edit','String',num2str(getfield(p,fnames{k})),'Position',pos + twid ,'BackGroundColor',[1 1 1], ...
                'Tag', ['mesoGT_pedit_' fnames{k}]);
             pos = pos + thei;         
         end;
     end;
     
     uicontrol('Style','pushbutton','String','OK','Position',[50 10 70 35],'Callback',{@OK,gcbo,[],[],[]});
     uicontrol('Style','pushbutton','String','Cancel','Position',[150 10 70 35],'Callback',{@CANCEL,gcbo,[],[],[]});

     uiwait(gcf);
     
     

function OK(h,eventdata,handles,varargin) 
    
    
    
     for k = 1:length(fnames),         
         if not(isa(getfield(p,fnames{k}),'function_handle')),
             h = findobj(gcf,'Tag', ['mesoGT_pedit_' fnames{k}]);
             val = str2num(get(h,'string'));
             p = setfield(p,fnames{k},val);
         end;
     end;  
   
    close(fignum);     
end

    
function CANCEL(h,eventdata,handles,varargin) 
    
    
    close(fignum);     
end

end
     
