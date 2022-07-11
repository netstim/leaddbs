function dicell=ea_dir2cell(di)
dicell={};
cnt=1;
for d=1:length(di)
    if ~strcmp(di(d).name(1),'.')
        dicell{cnt}=di(d).name;
        cnt=cnt+1;
    end
end