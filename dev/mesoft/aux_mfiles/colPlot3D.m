function colPlot3D(lines,col,csc)

N = 256;
cmap = autumn(N);
cmap = hot(N);
lw = 1;
delete(findobj(gca,'tag','fibplot3D_colplot3D'));
    for k = 1:length(lines),        
        fib = lines{k}';
        fib = fib';
        if size(fib,1) > 1,
            sp = 1;
            if nargin == 1,
                dif = (fib(2:end,:) - fib(1:end-1,:)).^2; dif = dif(:,[2 1 3]);
                dif = dif.*repmat(1./(eps+sum(dif,2)),[1 3]);        
                dif = [dif(1,:) ; dif];
                clear cold;
                cold(1,:,:) = dif; cold(2,:,:) = dif;
            else
                c = round(N*(col{k}-csc(1))/(csc(2)-csc(1)))+1;
                c(c<1) = 1; c(c>=N) = N;
                cold =repmat(reshape(cmap(c,:),[1 size(c,2) 3]),[2 1 1]);
            end;         
            handle = surface(ones(2,2),ones(2,2),ones(2,2),ones(2,2),'facecolor','no','edgecolor','interp','tag','fibplot3D_colplot3D','hittest','off','Clipping','off','linewidth',lw);
            set(handle,'xdata',[fib(1:sp:end,2) fib(1:sp:end,2)]','ydata',[fib(1:sp:end,1) fib(1:sp:end,1)]','zdata',[fib(1:sp:end,3) fib(1:sp:end,3)]','cdata',double(cold(:,1:sp:end,:)));
        end;
    end;