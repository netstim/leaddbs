function ea_show_SNR_res(allRes)


fn=fieldnames(allRes);

for sub=1:length(allRes)
    for f=1:length(fn)
        if isempty(allRes(sub).(fn{f}))
            allRes(sub).(fn{f}).c1snr=nan;
            allRes(sub).(fn{f}).c2snr=nan;
            allRes(sub).(fn{f}).c3snr=nan;
        end
    end

end

for f=1:length(fn)
    for sub=1:length(allRes)
        [c1vals(sub)]=allRes(sub).(fn{f}).c1snr;
        [c2vals(sub)]=allRes(sub).(fn{f}).c2snr;
        [c3vals(sub)]=allRes(sub).(fn{f}).c3snr;
    end

    h=figure;
    title(ea_underscore2space(fn{f}))
    ea_violin([c1vals',c2vals',c3vals']);
    g=gca;
    g.XTick=[1:3];
    g.XTickLabel={'SNR GM','SNR WM','SNR CSF'};
    axis square
    ylabel(['SNR ',ea_underscore2space(fn{f})]);
    saveas(h,['SNR_',(fn{f}),'.png']);
    close(h)
end
