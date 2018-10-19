function X=ea_graphvarmat2mat(txtfname,parcellation)

parcfid=fopen([ea_space([],'labeling'),parcellation,'.txt']);
labA=textscan(parcfid,'%f %s\n');
labels=labA{2};
X=nan(length(labels));


d=dir(txtfname);
for fi=1:length(d)
    fid=fopen(d(fi).name);
    A=textscan(fid,'%s');
    M=cellfun(@str2double,A{1});
    toprow=find(~isnan(M));
    Nx=toprow(1)-2;
    
    usedlabels=A{1}(1:Nx);
    [~,usedidx]=ismember(usedlabels,labels);
    
    
    
    M(1:Nx)=[];
    M(1:Nx+1:end)=[];
    M=reshape(M,length(usedlabels),length(usedlabels));
    
    X(usedidx,usedidx)=M;
end


