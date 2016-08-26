function ifcell=ea_checknoerrorfolders(ifcell)
cnt=1;
ecnt=1;

for f=1:length(ifcell)
   
    errfiles=dir([ifcell{f},filesep,'*.err']);
    
    for e=1:length(errfiles)
        fid = fopen([ifcell{f},filesep,errfiles(e).name]);
        A=textscan(fid,'%s');
        A=A{1};
        fclose(fid);
        if length(A)<2 % no error happened
            noerror(ecnt)=f;
            ecnt=ecnt+1;
        end
    end
    cnt=cnt+1;
end
try
ifcell(noerror)=[];
end