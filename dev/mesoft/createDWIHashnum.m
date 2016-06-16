function hashnum = createDWIHashnum(datastruct,Pstruc)


ten = datastruct.original_bTensor;
alpha = Pstruc.alpha + Pstruc.restrictions*2 + 8*Pstruc.sphericalDiscNumber;
ordermax = Pstruc.ordermax;
tissue = datastruct.fixedTissueMode;

for k = 1:5,
    hashnum(k) = sum(ten(:).^k);
end;

if isempty(tissue) 
    tissue = [-1 -1 -1];
end;

W = diag(createWeightingScheme(ten,Pstruc.b_weighting));
for k = 1:3,
    hashnum2 = sum(W.^k);
end;

hashnum = [alpha ordermax hashnum tissue hashnum2];
