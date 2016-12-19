function [res, errStr]= ten2camino(bTensor, scheme_fName)
%function [res, errStr]= ten2camino(bTensor, scheme_fName)
%
%
% Bjoern W. Kreher
% 04/08
%
% UNIX

res= []; errStr= '';

diffTime= 10^-1; % [s]


[bValAy, dirVc]= ten2dir(bTensor);

imNo= size(bTensor, 3);
waveVc= dirVc.*(sqrt(bValAy/diffTime)'*ones(1, 3))*1000;

[fHd, errStr]= fopen(scheme_fName, 'w+t');
if fHd < 0
    return;
end

fprintf(fHd, '%g\n%g\n', diffTime, imNo);

for i= 1:imNo
   fprintf(fHd, '%g\n%g\n%g\n', waveVc(i, 1), waveVc(i, 2), waveVc(i, 3)); 
end


fclose(fHd);
res= scheme_fName;




