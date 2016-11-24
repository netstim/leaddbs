function stackview(hr,showfun)



if nargin == 1,
    showfun = @(x) imagesc(x,[0 1]);
end;

[ds.maps cmap] =  compMesoParams(hr);
ds.img = ds.maps(:,:,:,4);
ds.hr = hr;
dirs = load('dwidirections');
ds.dirs256 = [dirs.dirs256 -dirs.dirs256];
ds.modelstruc = evalin('base','modelstruc');

proto.ten = hr.user.bTensor;
proto.b = squeeze(round((proto.ten(1,1,:)+proto.ten(2,2,:)+proto.ten(3,3,:))));     
proto.buni = unique(round(proto.b/100))*100;  
for k = 1:size(proto.ten,3);
    
    [U D] = eigs(proto.ten(:,:,k));
    [~,ix]=sort(D(logical(eye(length(D)))));
    proto.dirs(:,k) = U(:,ix(3));
    
end;

ds.proto = proto;






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

    
    b = ds.proto.b;
     
    Smeas = squeeze(ds.hr.dataAy(c(1),c(2),c(3),:));
    Smeas = Smeas ./ mean(Smeas(round(b/100)==0));

    buni = unique(round(b/100))*100;
  
        for k = 1:size(ds.hr.user.bTensor,3);
            
            [U D] = eigs(ds.hr.user.bTensor(:,:,k));
            [~,ix]=sort(D(logical(eye(length(D)))));
            dirs(:,k) = U(:,ix(3));

        end;

        
        %% get parameters of segment
                        
        Dpara_in = ds.maps(c(1),c(2),c(3),1);
        Dpara_ex = ds.maps(c(1),c(2),c(3),2)+ds.maps(c(1),c(2),c(3),3);
        Dorth =  ds.maps(c(1),c(2),c(3),3);
        vf =        ds.maps(c(1),c(2),c(3),4);
        vf(vf>1) = 1;
        vf_csf = ds.maps(c(1),c(2),c(3),5);
        vf_csf(vf_csf<0) = 0;
        snr = ds.maps(c(1),c(2),c(3),7);
        
        
        
        disptypes = {'poisson','watson','mises','heat'};
        res = fitDispersion(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,snr,ds.modelstruc.noisedeg, ds.proto, Smeas,ds.dirs256,disptypes);
        
           
        subplot(2,4,4);
        cla;
        n = res.n;
        for k = 1:size(n,2), 
            w = res.weights(k);
            plot3([n(1,k) -n(1,k)]*w,[n(2,k) -n(2,k)]*w,[n(3,k) -n(3,k)]*w,'linewidth',2); hold on;
        end;
        pd = res.pd;
       plot3([pd(1) -pd(1)],[pd(2) -pd(2)],[pd(3) -pd(3)],'r');
        
       showGlyph(ds.dirs256,res.fod*1);
        
        hold off;
        axis equal;
        axis([-1 1 -1 1 -1 1]);
        axis off;
        %view(-50,90);
             
%        subplot(2,4,4);
%         cla;    
%         
%         showGlyph(ds.dirs256,res.fod);
%         
%         hold off;
%         axis equal;
%         axis([-1 1 -1 1 -1 1]);
%         axis off;
        
 %%   
    subplot(2,4,[7 8]);

    style = {'-','-.','.','--'};
    
    col = 'rgbmc';
    colcnt = 1;
    for k = 2:length(buni),
        bidx = find(round(b/100) == round(buni(k)/100));
        
        proj = squeeze(cmult(res.pd,cmult(res.pd,ds.hr.user.bTensor(:,:,bidx),1,1),1,2))/buni(k);        
        [proj sidx] = sort(proj);
        
        for j = 1:length(disptypes),
            plot(proj,res.S{j}(bidx(sidx)),[style{j} col(colcnt)]); hold on;
        end;
        
        plot(proj,(squeeze(Smeas(bidx(sidx)))),['*' col(colcnt)]); hold on;

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
    


          subplot(2,4,3);
        cla;
        axis([0 10 0 10]);
        str = sprintf('Numcomps: %i \n',sum(res.weights>0));
        str = [str sprintf('vfintra : %.3f \n ',vf(1))];
        str = [str sprintf('vfcsf : %.3f  \n',vf_csf(1))];
        str = [str 'Dparallel:  ' sprintf('%.2f ',Dpara_in(1)) sprintf('\n')];
        str = [str 'Dparallel ex:  ' sprintf('%.2f ',Dpara_ex(1)) sprintf('\n')];
        str = [str 'Dortho:  ' sprintf('%.2f ',Dorth(1)) sprintf('\n')];
        str = [str 'dispersion kappa:  ' sprintf('%.2f ',res.kappa) sprintf('\n')];
        str = [str 'snr:  ' sprintf('%.2f ',snr(1)) sprintf('\n')];
        for k = 1:length(disptypes),
            str = [str disptypes{k} ' loglik (chi/gau):  ' sprintf('%.2f/%.2f ',res.errRic(k),res.errGS(k)) sprintf('\n')];
        end;
        text(1,5,str,'fontsize',25);
        axis off;
        
         output_txt = sprintf('xx');
%             
%      
% function [M idx] = SH(n,lmax)
% n = n ./ repmat(sqrt(sum(n.^2)),[3 1]);
% 
% M = [];
% for l = 0:2:lmax,
%     m = legendre(l,n(3,:),'sch');
%     n2 = n(1,:)+i*n(2,:);
%     n2 = n2 ./ abs(n2);
%     m = m.* (repmat(n2,[size(m,1) 1]).^repmat((0:l)',[1 size(n2,2)]))*sqrt(2*l+1);
%     idx1 = size(M,1);
%     M = [M ; m(1,:) ; real(m(2:end,:)) ; imag(m(2:end,:))];
%     idx2 = size(M,1);
%     idx{l/2+1} = idx1+1:idx2;
% end;
% 
% M = M/sqrt(size(n,2));
% 
% 
%  function S = simSig(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,ten,n,weights,disptype,kappa)
%  S = 0;
%  tr = squeeze(ten(1,1,:)+ten(2,2,:)+ten(3,3,:))/1000;
%  for k = 1:size(n,2),
%      q2 = squeeze(cmult(n(:,k),cmult(n(:,k),ten,1,1),1,2))/1000;
%      if strcmp(disptype,'nodisp') ~= 1,
%          Eintra = dispstick(Dpara_in,kappa,[],[q2(:) (tr(:)-q2(:))]',disptype);
%          Eextra = dispstick(Dpara_ex-Dorth,kappa,[],[q2(:) (tr(:)-q2(:))]',disptype) .* exp(-(tr-q2)*Dorth); 
%      else
%         Eintra = exp(-q2*Dpara_in);
%         Eextra = exp(-q2*Dpara_ex-(tr-q2)*Dorth); 
%      end;
%      Ecsf = exp(-tr*3);
%      S = S + weights(k)*(vf*Eintra+ (1-vf-vf_csf)*Eextra + vf_csf*Ecsf );
%  end;
% 
%  
%  
%  
%  
%  
%  
%  
% function res = fitDispersion(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,snr, ten,Smeas,dirs256,disptypes)
%  
%     
%     b = squeeze(round((ten(1,1,:)+ten(2,2,:)+ten(3,3,:))));
%      
% 
%     buni = unique(round(b/100))*100;
%   
%     for k = 1:size(ten,3);
%         [U D] = eigs(ten(:,:,k));
%         dirs(:,k) = U(:,1);
%     end;
%     
%     
%     
%     
%     
%     
%     
%        %% recompute signal
%         
%          lag = @(x) exp(x/2).*((1-x).*besseli(0,-x/2) - x.* besseli(1,-x/2));
%          riceexp = @(x,sigma) sigma*sqrt(pi/2) *lag(-x.^2./(2*sigma^2));        
% 
%       
%          S = simSig(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,ten,[0 0 1]',1,'nodisp');
%         
%          highb = round(b/100) == round(buni(end)/100);
%          
%          lmax = 8;
%          [ n weights fod] =  csd(Smeas(highb),dirs(:,highb),dirs256, S(highb),lmax);
%   
%         weights(weights<0.25) = 0;
%         weights = weights/sum(weights);
%         
%         subplot(2,4,4);
%         cla;
%         Mten = 0;
%         for k = 1:size(n,2), 
%             w = weights(k);
%             plot3([n(1,k) -n(1,k)]*w,[n(2,k) -n(2,k)]*w,[n(3,k) -n(3,k)]*w,'linewidth',2); hold on;
%             Mten = Mten + n(:,k)*n(:,k)' *w;
%         end;
%         [U D] = eigs(Mten);
%         pd = U(:,1);
%         plot3([pd(1) -pd(1)],[pd(2) -pd(2)],[pd(3) -pd(3)],'r');
%         
%         showGlyph(dirs256,fod*1);
%         
%         hold off;
%         axis equal;
%         axis([-1 1 -1 1 -1 1]);
%         axis off;
%         %view(-50,90);
%        
%          
%       %%       
%       
%       
%         lmax = 8;
%         kappa = 1;
%         S = simSig(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,ten,n,weights,'nodisp');
%         S = riceexp(S,1/snr);
%         for k = 2:length(buni),
%             bidx = round(b/100) == round(buni(k)/100);                  
%             [M shidx] = SH(dirs(:,bidx),lmax);
%             pM = M*S(bidx); 
%             pp = M*Smeas(bidx); 
%             for l = 2:length(shidx)
%                 pMP(l-1,k-1) = pp(shidx{l})'*pM(shidx{l});
%                 pPred2(l-1,k-1) = sum(pM(shidx{l}).^2);
%                 pMeas2(l-1,k-1) = sum(pp(shidx{l}).^2);
%             end;
%         end;
%                          
%         aa = sum(pMP,2)./sum(pPred2,2);
%         kappa = aa(1);
%         
%         
% 
%         
%         
%         
%         
%     clear S
%     for k = 1:length(disptypes);
%         S{k} = simSig(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,ten,n,weights,disptypes{k},kappa);  
%     end;
%     
%         
%      errRic = RicianLogLik(double(Smeas), double(S{1}),double( 1/snr(1)));
%      errGS = sum((Smeas-S{1}).^2)*snr(1)^2;
% 
%     for k = 1:length(disptypes);
%         S{k} = riceexp(S{k},1/snr(1));
%     end;
%  
%  
%     res.S = S;
%     res.errRic = errRic;
%     res.errGS = errGS;
%     res.fod = fod;
%     res.weights = weights;
%     res.n = n;
%     res.pd = pd;
%     res.kappa = kappa;
%     
 
 
 
 
 
 
 
 
 
 
            