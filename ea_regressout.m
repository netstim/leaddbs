function [Cleaned,beta_hat,nnanix]=ea_regressout(S,Cova)
% simple function that regresses out one variable from another



nnanix=any(isnan(Cova),1);
nnanix=nnanix+(isnan(S));
nnanix=~nnanix;
Cleaned=S;
Cova=Cova(:,nnanix);

S=S(nnanix);

Cova=zscore(Cova')';
S=zscore(S);

beta_hat = (Cova*Cova')\Cova*squeeze(S'); % ordinary least-squares estimator
Cleaned(nnanix)=squeeze(S)'-Cova'*beta_hat; % regress out each covariate.
Cleaned=Cleaned';
