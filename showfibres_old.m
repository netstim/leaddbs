function fib_plots=showfibres(fig,coords,options)



load(fullfile('fibres','gibbsconnectome10.mat'));
cnt=1;
for fib=1:100:length(gibbsconnectome)    
    
    for contact=1:8
        [IDX,D]=rangesearch(coords(contact,:),10,gibbsconnectome{fib},2);
        %plot3(gibbsconnectome{fib}(:,1),gibbsconnectome{fib}(:,2),gibbsconnectome{fib}(:,3),'k-');
        if ~isempty(IDX)
        connectingfibs{cnt}=gibbsconnectome{fib};
        cnt=cnt+1;
        end
    end
end
if ~isempty(connectingfibs)
    plot=1;
for fib=1:length(connectingfibs)
    for segment=1:length(connectingfibs{fib})-1;
    segclr=detcolor(connectingfibs{fib}(segment:segment+1,:));
    fib_plots(plot)=plot3(connectingfibs{fib}(segment:segment+1,1),connectingfibs{fib}(segment:segment+1,2),connectingfibs{fib}(segment:segment+1,3),'Color',segclr);
   plot=plot+1;
    end
    end
end

function rgb=detcolor(mat) % determine color based on traversing direction.

xyz=abs(diff(mat));
rgb=xyz/max(xyz(:));

