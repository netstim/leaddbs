
%%%% 
% specifiy here the norm used for the likelihood 
% chi2 = (S-M)'*W*(S-M)
%%%%

function [W Q] = createWeightingScheme(ten,bweight)

% be careful: not really included in the hashnum (sometimes there is an ambig.) so you have to rm -r /tmp/mesoGT/*

b = squeeze(ten(1,1,:) + ten(2,2,:) + ten(3,3,:));

% for the modelling term
b1 = b>0.1;
b0 = b<=0.1;


W = bweight * double(b0)/sum(b0) + double(b1)/sum(b1);

W = diag(W);
W = W / trace(W); % lot's of computation rely on this normalization, so keep it!!


% for the tracking guide
b1 = b>max(b)*0.9;
Q =  double(b1);
Q = diag(Q);

Q = Q /trace(Q);
