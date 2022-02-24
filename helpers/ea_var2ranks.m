function ranks=ea_var2ranks(var)
% function to convert variables to ranked data.


        [varsort,idx]=sort(var);
        ranks=zeros(length(varsort),1);
        for rank=1:length(ranks)
            ranks(idx(rank))=rank;
        end