function [select,list]=ea_groupselectorwholelist(select,list)


selpts=length(select);
allpts=length(list);
if ~(selpts==allpts)
    
    answ=questdlg(['You selected ',num2str(selpts),' from the whole group of ',num2str(allpts),' patients. Do you wish to perform the process on the selection or on all patients?'],'Process selection or whole set','Selection','Whole Group','Whole Group');
    
    switch lower(answ)
        case 'selection'
            return
        case 'whole group'
            select=1:length(list);
    end
end


