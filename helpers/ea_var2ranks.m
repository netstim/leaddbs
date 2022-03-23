function ranks=ea_var2ranks(var,order)
% function to convert variables to ranked data.
if ~exist('order','var')
    order='ascend';
end

        [varsort,idx]=sort(var,order);
        ranks=zeros(length(varsort),1);
        for rank=1:length(ranks)
            ranks(idx(rank))=rank;
        end