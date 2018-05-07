function stackview(mr,hr,ftr,showfun)



if nargin == 3,
    showfun = @imagesc;
end;

ds.img = mr;
ds.hr = hr;
ds.ftr = ftr;


set(datacursormode,'UpdateFcn',{@(x,y) labeldtips(x,y,ds)})
    


%f = figure('Position',[360,500,750,585]);    
f = gcf;
clf;

h = uicontrol('Style','slider','String','Surf','Position',[0,0,300,25],'Tag','theslider','Value',0.5,'Callback',{@resolutionslider_Callback,gcbo,[],[],ds,showfun});
botPanel = uibuttongroup('Units','pixels','Position',[400 0  200 30],'SelectionChangeFcn',{@radiob_Callback,gcbo,ds},'Tag','butgroup');
hX = uicontrol('Style','radiobutton','String','X','Position',[10,1,40,25],'Tag','1','Parent',botPanel);
hY = uicontrol('Style','radiobutton','String','Y','Position',[60,1,40,25],'Tag','2','Parent',botPanel);
hZ = uicontrol('Style','radiobutton','String','Z','Position',[110,1,40,25],'Tag','3','Parent',botPanel);
coord = uicontrol('Style','text','String','x:1','Tag','coordinfo','Position',[300,0,100,25]);
set(h,'SliderStep',[1/size(ds.img,1) 0.5]);
showslice(ds.img,showfun);




% --------------------------------------------------------------------
function varargout = resolutionslider_Callback(h, eventdata, handles, varargin)
  
    ds = varargin{3};
    showslice(ds.img,varargin{4});
    
% --------------------------------------------------------------------
function varargout = radiob_Callback(h, eventdata, handles, varargin)
    img = varargin{1}.img;
     
    radiob = get(h,'SelectedObject');
    selection = str2num(get(radiob,'Tag'));    
    
    slider = findobj(gcf,'Tag','theslider');
    set(slider,'SliderStep',[1/size(img,selection) 0.5]);
    cb = get(slider,'Callback');
    
    showslice(img,cb{6});
    
function showslice(img,showfun);
    subplot(2,4,[1 2 5 6]);
    slider = findobj(gcf,'Tag','theslider');
    pos = get(slider,'value');
    botpanel = findobj(gcf,'Tag','butgroup');
    radiob = get(botpanel,'SelectedObject');
    selection = str2num(get(radiob,'Tag'));    
    idxpos = round(pos*(size(img,selection)-1))+1;
    cinfo = findobj(gcf,'Tag','coordinfo');
    
    switch selection,
        case 3,
            showfun(squeeze(img(:,:,idxpos,:)));
            set(cinfo,'String',sprintf('%i/%i',idxpos,size(img,3)));    
        case 2,
            showfun(flipdim(permute(squeeze(img(:,idxpos,:,:)),[2 1 3]),1));
            set(cinfo,'String',sprintf('%i/%i',idxpos,size(img,2)));    
        case 1,
            showfun(flipdim(permute(squeeze(img(idxpos,:,:,:)),[2 1 3]),1)); 
            set(cinfo,'String',sprintf('%i/%i',idxpos,size(img,1)));    
    end;
          
            
     
function output_txt = labeldtips(obj,event_obj,ds)
img = ds.img;
    slider = findobj(gcf,'Tag','theslider');
    pos = get(slider,'value');
    botpanel = findobj(gcf,'Tag','butgroup');
    radiob = get(botpanel,'SelectedObject');
    selection = str2num(get(radiob,'Tag'));    
    idxpos = round(pos*(size(img,selection)-1))+1;
    cinfo = findobj(gcf,'Tag','coordinfo');
    cursorpos = get(event_obj,'Position');    
    
    switch selection,
        case 3,
            c = [cursorpos([2 1]) idxpos];
        case 2,
            c = [cursorpos(1) idxpos size(img,3)-cursorpos(2)+1 ];
        case 1,
            c = [idxpos cursorpos(1) size(img,3)-cursorpos(2)+1];
    end;

    imval = ds.img(c(1),c(2),c(3));
    
    
    ftr = ds.ftr;
    osamp = (ftr.trackParam.params.p_wid);
     wscale = (ftr.trackParam.params.p_weight);
    pos = ftr.user.P(1:3,:)';
    pos = (pos ./ repmat(ftr.vox(:)',[size(pos,1) 1]));
    idx = find(floor(pos(:,1))+1 == c(1) & floor(pos(:,2))+1 == c(2) & floor(pos(:,3))+1 == c(3));
    
    b = round((ds.hr.user.bTensor(1,1,:)+ds.hr.user.bTensor(2,2,:)+ds.hr.user.bTensor(3,3,:))/100);
     
    Smeas = squeeze(ds.hr.dataAy(c(1),c(2),c(3),:));
    Smeas = Smeas ./ mean(Smeas(b==0));

    
    if not(isempty(idx)),
        
        pos = ftr.user.P(1:3,idx)'./ repmat(ftr.vox(:)',[length(idx) 1]);
        pos = floor(pos*osamp);
        idx_single = find(all(pos == repmat(pos(1,:),[size(pos,1) 1]),2));
        idx = idx(idx_single);
        pos = pos(idx_single,:);

        
        %% get parameters of segment
        n = ftr.user.P(4:6,idx);
%         
%         weights = ftr.user.P(14,idx)*wscale;
%         Dpara = ftr.user.P(11,idx);
%         Dorth = ftr.user.P(13,idx);
%         for k = 1:length(idx),
%             vfi(k) = ftr.user.vf(pos(k,1)+1,pos(k,2)+1,pos(k,3)+1,1);
%             psw(k) = ftr.user.vf(pos(k,1)+1,pos(k,2)+1,pos(k,3)+1,2);
%         end;
%         
%         %% recompute signal
%         W = createWeightingScheme(ds.hr.user.bTensor);
%         S = 0;
%         tr = squeeze(ds.hr.user.bTensor(1,1,:)+ds.hr.user.bTensor(2,2,:)+ds.hr.user.bTensor(3,3,:))/1000;
%         for k = 1:length(idx),
%             q2 = squeeze(cmult(n(:,k),cmult(n(:,k),ds.hr.user.bTensor,1,1),1,2))/1000;
%             Eintra = exp(-q2*Dpara(k));% Eintra = Eintra - sum(W*Eintra);
%             Eextra = exp(-q2*Dpara(k)-(tr-q2)*Dorth(k)); % Eextra = Eextra - sum(W*Eextra);
%             S = S + weights(k)*(vfi(k)*Eintra+ (1-vfi(k))*Eextra);
%         end;
%         mean_measS = ftr.user.meansignal(ftr.user.S2(c(1),c(2),c(3),2));
%         vfsw = (mean_measS-psw(1));
%         S = S + vfsw;
     
      %%
        
      
      
        
        subplot(2,4,4);
        Mten = 0;
        for k = 1:length(idx),    
            w = 1; %weights(k);
            plot3([n(1,k) -n(1,k)]*w,[n(2,k) -n(2,k)]*w,[n(3,k) -n(3,k)]*w,'linewidth',2); hold on;
            Mten = Mten + n(:,k)*n(:,k)' *w;
        end;
        [U D] = eigs(Mten);
        pd = U(:,1);
        plot3([pd(1) -pd(1)],[pd(2) -pd(2)],[pd(3) -pd(3)],'r');
        
        hold off;
        axis equal;
        axis([-1 1 -1 1 -1 1]);
        axis off;
%             
%       
%           subplot(2,4,3);
%         cla;
%         axis([0 10 0 10]);
%         str = sprintf('Numcomps: %i \n',length(Dpara));
%         vfi_abs = (weights.*vfi)./(vfsw+sum(weights));
%         str = [str sprintf('vfintra : %.3f  ',sum(vfi_abs)) '(' sprintf('%.2f ',vfi_abs) sprintf(')\n')  ];
%         vfe_abs = (weights.*(1-vfi))./(vfsw+sum(weights));
%         str = [str sprintf('vfextra : %.3f  ',sum(vfe_abs)) '(' sprintf('%.2f ',vfe_abs) sprintf(')\n')  ];
%         str = [str sprintf('vfstill : %.3f  ',vfsw)  sprintf('\n')  ];
%         str = [str 'weights : ' sprintf('%.3f ',weights)  sprintf('\n')  ];
%         str = [str sprintf('S0 : %.3f  ',vfsw+sum(weights))  sprintf('\n')  ];
%         str = [str 'Dparallel:  ' sprintf('%.2f ',Dpara) sprintf('\n')];
%         str = [str 'Dortho:  ' sprintf('%.2f ',Dorth) sprintf('\n')];
%         str = [str 'xi2:  ' sprintf('%.3f ',std(Smeas-S)) sprintf('\n')];
%         text(1,5,str,'fontsize',15);
%         axis off;
        
        
      
        output_txt = sprintf(['map val:' num2str(imval)]);
    else
        output_txt = sprintf('empty');    
    end;

    
    subplot(2,4,[7 8]);
    buni = unique(b);
    
    col = 'rgbmc';
    colcnt = 1;
    for k = 2:length(buni),
        bidx = b == buni(k);
        if exist('S'),
            Sbidx = squeeze(S(bidx));
            [Sbidx sidx] = sort(Sbidx);
            proj = squeeze(cmult(pd,cmult(pd,ds.hr.user.bTensor(:,:,bidx),1,1),1,2))/buni(k)/100;        
            plot(proj(sidx),Sbidx,['-' col(colcnt)]); hold on;
        else
           % pd = [1 0 0]';
            proj = squeeze(cmult(pd,cmult(pd,ds.hr.user.bTensor(:,:,bidx),1,1),1,2))/buni(k)/100;        
        end;        
        plot(proj,squeeze(Smeas(bidx)),['*' col(colcnt)]); hold on;
        colcnt = colcnt + 1;
    end;
    line([0 1],[0 0],'color','k');
    hold off;
    axis([0 1 -0.2 1.2]);
    grid on;
%    legend('b=1000 fit','b=1000 signal','b=2000 fit','b=2000 signal','b=3000 fit','b=3000 signal')
    legend('b=1000 fit','b=1000 signal','b=2000 fit','b=2000 signal')
     xlabel('cos(theta) to prinipal direction')
     %  output_txt = sprintf(['hallo ' num2str(imval)]);
    
    
            
            