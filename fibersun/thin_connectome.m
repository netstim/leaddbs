function thin_connectome


load('gibbsconnectome');
gc=gibbsconnectome;
for thin=[2]
    disp(['Thinning by ',num2str(thin),'.']);
    gibbsconnectome=gc(1:thin:end);
    save(['gibbsconnectome',num2str(thin)],'gibbsconnectome');
    
end