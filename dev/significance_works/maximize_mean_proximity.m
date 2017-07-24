load('XYZV');
for save=1:500
for iter=1:2
    
    dat=XYZV;
    if iter<2
        dat(:,4)=dat(randperm(length(dat)),4);
    end

Peuc=squareform(ea_pdist(dat(:,1:3))); % distance matrix
Pval=squareform(ea_pdist(dat(:,4))); % distance matrix

% exp
Peuc(logical(eye(size(Peuc,1))))=1./exp(Peuc(logical(eye(size(Peuc,1)))));
Pval(logical(eye(size(Peuc,1))))=1./exp(Pval(logical(eye(size(Peuc,1)))));


[X,Y]=meshgrid(dat(:,4),dat(:,4));
vals=mean(cat(3,X,Y),3);
Pval=vals.^(Pval);

P=Peuc.^Pval;
P(logical(eye(size(P,1))))=0;
P=sum(P);
P=P;
m(iter)=max(P);
mm(iter,save)=max(P);
end
di(save)=m(2)-m(1);
mi(save)=nanmean(m);
end
di=mean(di);
mi=mean(mi);
disp(['Distance difference: ', num2str((di/mi)*100),'%.']);
disp(['Rough std is: ',num2str(std(P(:))/mean(P(:))),'.']);
mean(mm')
