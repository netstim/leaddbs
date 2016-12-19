function result= calc_i(D, bTensor)
%function result= calc_i(D, bTensor)
%
%
% Bjoern W. Kreher
% 01/03
%
% UNIX

if ischar(bTensor)
    tmp= open(bTensor);
    b_tensor= tmp.B_tensor;
    clear tmp;
else
    b_tensor= bTensor;
end
gradNo= size(b_tensor, 3);

result= [];

for gradI= 1:1:gradNo
    result(gradI)= exp(-sum(reshape(D.*b_tensor(:, :, gradI), [numel(D) 1])));
end
