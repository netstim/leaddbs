function stackview(hr,ftr,showfun)



if nargin == 2,
    showfun = @imagesc;
end;

ds.maps = computeParameterMaps(ftr,1);
ds.img = ds.maps.vf;
ds.hr = hr;
ds.ftr = ftr;

set(datacursormode,'UpdateFcn',{@(x,y) labeldtips(x,y,ds)})
    
datacursormode on

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
        
        weights = ftr.user.P(14,idx)*wscale;
        weights = weights/sum(weights);
        for k = 1:length(idx),
            Dpara_in(k) = ds.maps.Din(pos(k,1)+1,pos(k,2)+1,pos(k,3)+1,1);
            Dpara_ex(k) = ds.maps.Dexax(pos(k,1)+1,pos(k,2)+1,pos(k,3)+1,1);
            Dorth(k) = ds.maps.Dexrad(pos(k,1)+1,pos(k,2)+1,pos(k,3)+1,1);
            vf(k) =        ds.maps.vf(pos(k,1)+1,pos(k,2)+1,pos(k,3)+1,1);
            vf_csf(k) = ds.maps.vf_csf(pos(k,1)+1,pos(k,2)+1,pos(k,3)+1,1);
            snr(k) = ds.maps.snr(pos(k,1)+1,pos(k,2)+1,pos(k,3)+1,1);
        end;
    
        
      
      
        
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
            
      
      
        
    

        
        
        
        %% recompute signal
        

         lag = @(x) exp(x/2).*((1-x).*besseli(0,-x/2) - x.* besseli(1,-x/2));
         riceexp = @(x,sigma) sigma*sqrt(pi/2) *lag(-x.^2./(2*sigma^2));        
          snr = snr*sqrt(0.5);
         
         S = 0;
         tr = squeeze(ds.hr.user.bTensor(1,1,:)+ds.hr.user.bTensor(2,2,:)+ds.hr.user.bTensor(3,3,:))/1000;
         for k = 1:length(idx),
             q2 = squeeze(cmult(n(:,k),cmult(n(:,k),ds.hr.user.bTensor,1,1),1,2))/1000;
             
             Eintra = exp(-q2*Dpara_in(k));
             Eextra = exp(-q2*Dpara_ex(k)-(tr-q2)*Dorth(k)); 
             Ecsf = exp(-tr*3);
             S = S + weights(k)*(vf(k)*Eintra+ (1-vf(k)-vf_csf(k))*Eextra + vf_csf(k)*Ecsf );
         end;
         
        S = riceexp(S,1/snr(1));
         
         
      %%        
        buni = unique(b);
          
      
        for k = 1:size(ds.hr.user.bTensor,3);
            [U D] = eigs(ds.hr.user.bTensor(:,:,k));
            dirs(:,k) = U(:,1);
        end;
        lmax = 2;
        p1 = 0;
        p2 = 0;
        for k = 2:length(buni),
            bidx = b == buni(k);                  
            [M ] = SH(dirs(:,bidx),lmax);
            p = M*S(bidx); p1 = p1+sum(p(2:end).^2);
            p = M*Smeas(bidx); p2 =p2+sum(p(2:end).^2);            
            SHt{k-1} = M;
            %SHidx{k-1} = idx;
        end;
        kappa = sqrt(log(p1/p2))/8;
        
        
        
        
        
        
        %% recompute signal
         S = 0;
         tr = squeeze(ds.hr.user.bTensor(1,1,:)+ds.hr.user.bTensor(2,2,:)+ds.hr.user.bTensor(3,3,:))/1000;
         for k = 1:length(idx),
             q2 = squeeze(cmult(n(:,k),cmult(n(:,k),ds.hr.user.bTensor,1,1),1,2))/1000;
             Eintra = dispstick(Dpara_in(k),kappa,[],[q2(:) (tr(:)-q2(:))]');
             Eextra = exp(-q2*Dpara_ex(k)-(tr-q2)*Dorth(k)); % Eextra = Eextra - sum(W*Eextra);
             Ecsf = exp(-tr*3);
             S = S + weights(k)*(vf(k)*Eintra+ (1-vf(k)-vf_csf(k))*Eextra + vf_csf(k)*Ecsf );
         end;        
         S = riceexp(S,1/snr(1));

        S = abs(real(S));
        
%         SHT = zeros(size(SHt{1},1)*(length(buni)-1), size(dirs,2));       
%         SHordidx = zeros(size(SHt{1},1)*(length(buni)-1),1);
%         for k = 2:length(buni),
%             bidx = b == buni(k);                  
%             SHT((1:size(SHt{1},1))+(k-2)*size(SHt{1},1),bidx) = SHt{k-1};
%             for l = 0:2:lmax,
%                 SHordidx(SHidx{k-1}{l/2+1}+(k-2)*size(SHt{1},1)) = l/2+1;
%             end
%             
%         end;
%          l = (0:2:lmax)';
%          sm = exp(-0.005*l.*(l+1));
%          sm
%           d1 = SHT*S;
%           d2= sm(SHordidx).*d1;
%           S = pinv(SHT)*d2;
%         
%         sm(SHordidx)

%%





%%        
        

          subplot(2,4,3);
        cla;
        axis([0 10 0 10]);
        str = sprintf('Numcomps: %i \n',length(Dpara_in));
        str = [str sprintf('vfintra : %.3f \n ',vf(1))];
        str = [str sprintf('vfcsf : %.3f  \n',vf_csf(1))];
        str = [str 'Dparallel:  ' sprintf('%.2f ',Dpara_in(1)) sprintf('\n')];
        str = [str 'Dparallel ex:  ' sprintf('%.2f ',Dpara_ex(1)) sprintf('\n')];
        str = [str 'Dortho:  ' sprintf('%.2f ',Dorth(1)) sprintf('\n')];
        str = [str 'dispersion kappa:  ' sprintf('%.2f ',kappa) sprintf('\n')];
        str = [str 'snr:  ' sprintf('%.2f ',snr(1)) sprintf('\n')];
        text(1,5,str,'fontsize',15);
        axis off;
        
        
      
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
            plot(proj(sidx),(Sbidx),['-' col(colcnt)]); hold on;
        else
           % pd = [1 0 0]';
            proj = squeeze(cmult(pd,cmult(pd,ds.hr.user.bTensor(:,:,bidx),1,1),1,2))/buni(k)/100;        
        end;        
        plot(proj,(squeeze(Smeas(bidx))),['*' col(colcnt)]); hold on;
        colcnt = colcnt + 1;
    end;
    line([0 1],[0 0],'color','k');
    hold off;
    axis([0 1 0 1]);
    grid on;
%    legend('b=1000 fit','b=1000 signal','b=2000 fit','b=2000 signal','b=3000 fit','b=3000 signal')
    legend('b=1000 fit','b=1000 signal','b=2000 fit','b=2000 signal')
     xlabel('cos(theta) to prinipal direction')
     %  output_txt = sprintf(['hallo ' num2str(imval)]);
    
    
function [M idx] = SH(n,lmax)
n = n ./ repmat(sqrt(sum(n.^2)),[3 1]);

M = [];
for l = 0:2:lmax,
    m = legendre(l,n(3,:),'sch');
    n2 = n(1,:)+i*n(2,:);
    n2 = n2 ./ abs(n2);
    m = m.* (repmat(n2,[size(m,1) 1]).^repmat((0:l)',[1 size(n2,2)]));
    idx1 = size(M,1);
    M = [M ; m(1,:) ; real(m(2:end,:)) ; imag(m(2:end,:))];
    idx2 = size(M,1);
    idx{l/2+1} = idx1+1:idx2;
end;

M = M/sqrt(size(n,2));

            
            