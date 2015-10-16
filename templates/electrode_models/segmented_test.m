clear


h=figure;
hold on
for conins=1:2
    de=figure('visible','off');
    hold on
    switch conins
        case 1 % contacts are plotted
            ix=[15,105
                135,225
                255,345];
            fc='r';
        case 2 % insulations are plotted
            ix=[1,15
                105,135
                225,255
                345,360];
            fc='b';
    end
    
    
    for out=[1,2]
        switch out
            case 1 % small circle inside electrode
                r=0.5;
            case 2 % outside electrode
                r=1;
        end
        [X,Y,Z]=cylinder(r,360); % starts at 3 o clock and goes counter-clockwise.
        
        for i=1:size(ix,1)
            
            if conins==2 && out==1 % small cylinder
                smcy=surf(X,Y,Z);
                smcy=surf2patch(smcy,'triangles');
                smcy.vertices=[0,0,0
                    0,0,1
                    smcy.vertices];
                for smci=3:length(smcy.vertices)-2
                    smcy.faces=[smcy.faces;
                        smci,smci+2,1+(~mod(smci,2)*1)];
                end
            end
            
            
            
            segX{i}{out}=X(:,ix(i,1):ix(i,2));
            segY{i}{out}=Y(:,ix(i,1):ix(i,2));
            segZ{i}{out}=Z(:,ix(i,1):ix(i,2));
            oute{i}{out}=surf(segX{i}{out},segY{i}{out},segZ{i}{out});
            oute{i}{out}=surf(segX{i}{out},segY{i}{out},segZ{i}{out});
            oute{i}{out}=surf(segX{i}{out},segY{i}{out},segZ{i}{out});
            fv{i}{out}=surf2patch(oute{i}{out},'triangles');
        end
    end
    close(de)
    
    
    for i=1:size(ix,1)
        fvc{i}.faces=[fv{i}{1}.faces;fv{i}{2}.faces+length(fv{i}{1}.vertices)];
        
        fvc{i}.vertices=[fv{i}{1}.vertices;fv{i}{2}.vertices];
        
        % add plates..
        
        % figure
        % patch(fvc{i})
        fvc{i}.faces=[fvc{i}.faces
            2,1,length(fv{i}{1}.vertices)+2
            1,length(fv{i}{1}.vertices)+1,length(fv{i}{1}.vertices)+2];
        fvc{i}.faces=[fvc{i}.faces
            length(fv{i}{1}.vertices),length(fv{i}{1}.vertices)-1,length(fvc{i}.vertices)
            length(fvc{i}.vertices),length(fvc{i}.vertices)-1,length(fv{i}{1}.vertices)-1];
        
        
        
        
        %close lids (top/bottom)..
        for face=1:length(fv{i}{1}.vertices)-2;
            fvc{i}.faces=[fvc{i}.faces;
                face+length(fv{i}{1}.vertices)+2,face+2,face;
                face+length(fv{i}{1}.vertices)+2,face+length(fv{i}{1}.vertices),face];
        end
        
        set(0,'CurrentFigure',h);
        p{conins,i}=patch(fvc{i},'FaceColor',fc,'EdgeColor','none','CDataMapping','direct','CData',repmat([0.5,0.5,0.5],length(fvc{i}.faces),1));
        
        if conins==2 % plot ring of insulation
            
            p{conins,i+1}=patch(smcy,'FaceColor',fc,'EdgeColor','none','CDataMapping','direct','CData',repmat([0.5,0.5,0.5],length(smcy.faces),1));
        end
        
    end
    clear segX segY segZ oute fv fvc
end

a=camlight('headlight');
