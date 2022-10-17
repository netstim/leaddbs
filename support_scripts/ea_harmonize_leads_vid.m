% Example-script to create a video that iteratively migrates electrodes of a lead
% group scene to a common target.
% (c) 2022 Andreas Horn, BWH / HMS

daObj=VideoWriter('lg_video.mp4','MPEG-4'); % preferred profile

resultfig = gcf;

h=getappdata(resultfig,'el_render');

lcnt=1; rcnt=1;

for el=1:length(h)

    thishead=h(el).elstruct.markers(h(el).side).head;
    thistail=h(el).elstruct.markers(h(el).side).tail;
    thisx=h(el).elstruct.markers(h(el).side).x;
    thisy=h(el).elstruct.markers(h(el).side).y;

    switch h(el).side
        case 1
            heads.r(rcnt,:)=thishead;
            tails.r(rcnt,:)=thistail;
            x.r(rcnt,:)=thisx;
            y.r(rcnt,:)=thisy;
            rcnt=rcnt+1;
        case 2
            heads.l(lcnt,:)=thishead;
            tails.l(lcnt,:)=thistail;
            x.l(lcnt,:)=thisx;
            y.l(lcnt,:)=thisy;
            lcnt=lcnt+1;
    end
end

avghead.l=mean(heads.l);
avghead.r=mean(heads.r);
avgtail.l=mean(tails.l);
avgtail.r=mean(tails.r);
avgx.l=mean(x.l);
avgx.r=mean(x.r);
avgy.l=mean(y.l);
avgy.r=mean(y.r);


% compute vectors to apply in each iter
lcnt=1; rcnt=1;
for el=1:length(h)
    thishead=h(el).elstruct.markers(h(el).side).head;
    thistail=h(el).elstruct.markers(h(el).side).tail;
    thisx=h(el).elstruct.markers(h(el).side).x;
    thisy=h(el).elstruct.markers(h(el).side).y;

    switch h(el).side
        case 1
            vec.head.r(rcnt,:)=avghead.r-thishead;
            vec.tail.r(rcnt,:)=avgtail.r-thistail;
            vec.x.r(rcnt,:)=avgx.r-thisx;
            vec.y.r(rcnt,:)=avgy.r-thisy;
            rcnt=rcnt+1;
        case 2
            vec.head.l(lcnt,:)=avghead.l-thishead;
            vec.tail.l(lcnt,:)=avgtail.l-thistail;
            vec.x.l(lcnt,:)=avgx.l-thisx;
            vec.y.l(lcnt,:)=avgy.l-thisy;
            lcnt=lcnt+1;
    end
end

% apply vectors to electrodes:
Nsteps=100;
steps=exp(-0.05*[1:Nsteps]);
steps=steps/sum(steps); % together 1 step
open(daObj);
for initframes=1:100
        writeVideo(daObj,getframe(resultfig)); %use figure, since axis changes size based on view
end
ea_dispercent(0,'Writing video');
for iter=1:Nsteps
    writeVideo(daObj,getframe(resultfig)); %use figure, since axis changes size based on view
    lcnt=1; rcnt=1;
    for el=1:length(h)
        switch h(el).side
            case 1
                h(el).elstruct.markers(h(el).side).head=h(el).elstruct.markers(h(el).side).head+...
                    vec.head.r(rcnt,:)*steps(iter);
                h(el).elstruct.markers(h(el).side).tail=h(el).elstruct.markers(h(el).side).tail+...
                    vec.tail.r(rcnt,:)*steps(iter);
                h(el).elstruct.markers(h(el).side).x=h(el).elstruct.markers(h(el).side).x+...
                    vec.x.r(rcnt,:)*steps(iter);
                h(el).elstruct.markers(h(el).side).y=h(el).elstruct.markers(h(el).side).y+...
                    vec.y.r(rcnt,:)*steps(iter);
                rcnt=rcnt+1;
            case 2
                h(el).elstruct.markers(h(el).side).head=h(el).elstruct.markers(h(el).side).head+...
                    vec.head.l(lcnt,:)*steps(iter);
                h(el).elstruct.markers(h(el).side).tail=h(el).elstruct.markers(h(el).side).tail+...
                    vec.tail.l(lcnt,:)*steps(iter);
                h(el).elstruct.markers(h(el).side).x=h(el).elstruct.markers(h(el).side).x+...
                    vec.x.l(lcnt,:)*steps(iter);
                h(el).elstruct.markers(h(el).side).y=h(el).elstruct.markers(h(el).side).y+...
                    vec.y.l(lcnt,:)*steps(iter);
                lcnt=lcnt+1;
        end
    end
    drawnow
    ea_dispercent(iter/Nsteps);
end
    writeVideo(daObj,getframe(resultfig)); %use figure, since axis changes size based on view

ea_dispercent(1,'end');

close(daObj);

