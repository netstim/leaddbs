function [b_valErg, dirVc, b0Idx, bValIdx]= ten2dir(bTensor)
%function [b_val, dirVc, b0Idx, bValIdx]= ten2dir(bTensor)
%
%
% Bjoern W. Kreher
% 03/03
%
% UNIX


sizeAy= size(bTensor);

base=   [1 1 1]'*(0:9:(sizeAy(3)*9 - 1));

offset= [1 5 9]'*ones(1, sizeAy(3)); % [1 5 9] <=> [(1, 1) (2, 2) (3, 3)]
diagIdx= reshape(base + offset, [1 sizeAy(3)*3]);
diagElem= reshape(bTensor(diagIdx), [3 sizeAy(3)]);
b_valErg= sum(diagElem, 1);
b_val= b_valErg;
b_val= (b_val == 0) + b_val;

offset= [6 3 2]'*ones(1, sizeAy(3)); % [6 3 2] <=> [(1, 3) (2, 3) (1, 2)]
crossIdx= reshape(base + offset, [1 sizeAy(3)*3]);
crossElem= reshape(bTensor(crossIdx), [3 sizeAy(3)])./reshape(ones(3, 1)*b_val, [3 sizeAy(3)]);

sig= sign(crossElem);
idx= sum(abs(sig)) == 0;      %%% Fall 1: zwei Koponenten sind null
sig(:, idx)= 1;
modIdx= [1:3 1:3];                  %%% Fall 2: eine koponente ist null
idx= sum(abs(sig)) == 1;
sig(:, idx)= abs(abs(sig(:, idx)) - 1);
[dummy, idx2]= max(sig(:, idx));
if ~isempty(idx2)
    sig((idx - 1)*3 + idx2)= sign(sum(crossElem(:, idx)));
end

dirVc= (sqrt(diagElem./reshape(ones(3, 1)*b_val, [3 sizeAy(3)])).*sig)';

% dirVc= zeros(sizeAy(3), 3);
% 
% [dummy, idx]= max(
% dirVc(:, 3)= sqrt(M(2, 3, :).*M(1, 3, :)./M(1, 2, :));
% dirVc(:, 2)= M(2, 3, :)./dirVc(:, 3);
% dirVc(:, 1)= M(1, 2, :)./dirVc(:, 2);
% 
%dirVc= sqrt(M');

b0Idx= find(b_val < 100);
bValIdx= find(b_val >= 100);
