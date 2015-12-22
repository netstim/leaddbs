function smODF = smoothODF(ODF,bDir,numit)


smODF = ODF;
for k = 1:size(ODF,1);
    D = eye(3)*0.2 + single(bDir(:,k)*bDir(:,k)'); 
    smODF(k,:,:,:) = anisoDiffusionHomogenous(single(squeeze(ODF(k,:,:,:))),[D(1,1) D(2,2) D(3,3) D(1,2) D(1,3) D(2,3)]',single([0.1 numit]));
    fprintf('.');
end;
fprintf('\n');