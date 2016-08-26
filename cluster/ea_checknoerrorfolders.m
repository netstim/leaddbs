function ifcell=ea_checkforerrors(ifcell)


for f=1:length(ifcell)
   
    errfiles=dir([ifcell{f},'*.err']);
    for e=1:length(errfiles)
        fid = fopen([ifcell{f},errfiles(e).name]);
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
ifcell(noerror)=[];