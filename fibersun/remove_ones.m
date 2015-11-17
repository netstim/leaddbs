function remove_ones

for thin=[2,5,10,50,100,500]
        disp(['Thin is ',num2str(thin),'.']);

    load(['gibbsconnectome',num2str(thin)]);
    for fib=1:length(gibbsconnectome)
       gibbsconnectome{fib}=gibbsconnectome{fib}(:,1:3);
    end
    save(['gibbsconnectome',num2str(thin)],'gibbsconnectome');

end

load('gibbsconnectome');
    for fib=1:length(gibbsconnectome)
       gibbsconnectome{fib}=gibbsconnectome{fib}(:,1:3);
    end
    save(['gibbsconnectome'],'gibbsconnectome','-v7.3');
    