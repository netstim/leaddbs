function showTractStats(datastruct)

          ftr = datastruct.ftr;
          if isempty(ftr);
              return;
          end;
          fh = findobj('tag','mesoGT_statfigure');                            
          ssds = get(fh,'userdata');
          sfigure(fh);
          if isempty(findobj(gcf,'Tag','mesoGT_slicepopup')),                 
              uicontrol('Style','popupmenu','String', {'X','Y', 'Z'},'units','normalized','Position',[0.1 0.88 0.03 0.05],'Tag','mesoGT_slicepopup', 'Value',3,'Callback',{@combobox_Callback,gcbo,[],[]});       
              uicontrol('Style','slider','units','normalized','Position',[0.1 0.77 0.03 0.13],'Tag','mesoGT_slider','Value',0.5,'Callback',{@slider_Callback,gcbo,[],[]});                  
              uicontrol('Style','pushbutton','units','normalized','String','Browse Signal','Position',[0.03 0.03 0.2 0.05],'Tag','mesoGT_browsesig','Callback',{@browsesig_Callback,gcbo,[],[]});
              
          end;
          
          showFibs(ftr,ssds.fibcoloring);
          
           maps = computeParameterMaps(ftr,1);                            
          sfigure(fh);
          ssds.maps = cat(4,maps.Din,maps.Dexax,maps.Dexrad,maps.vf,maps.vf_csf,maps.segcount, maps.termcount,maps.snr);
          
          ssds.ds = datastruct;
          
          plotHists(datastruct,ssds.maps)
          
          
          set(fh,'userdata',ssds);
          showMaps;
          combobox_Callback([],[],[]) ;
          drawnow;
              
    
function showFibs(ftr,coloring)
         fh = findobj('tag','mesoGT_statfigure');                            
         sfigure(fh);
         subplot(4,4,[3 4 7 8]);
          numfibs = 1000;
          n = size(ftr.curveSegCell,1);
          if n>0,
              idx = rand(n,1)>(1-numfibs/n);
              if sum(idx)>0,                  
                  if coloring == 1,
                      colPlot3D(ftr.curveSegCell(idx));      
                  elseif coloring <= 4
                      colPlot3D(ftr.curveSegCell(idx),cellfun(@(x) x(coloring-1,:),ftr.curveD(idx),'uniformoutput',false),[0 3]); 
                  else
                      colPlot3D(ftr.curveSegCell(idx),cellfun(@(x) x(coloring-1,:),ftr.curveD(idx),'uniformoutput',false),[0 1]); 
                  end;                  
                  grid on;
                  set(gca,'PlotBoxAspectRatio',[1 1 1]);
                  set(gca,'Xlim',[-inf inf]);
                  set(gca,'Ylim',[-inf inf]);
                  set(gca,'Zlim',[-inf inf]);
              end;
          end;    
          contrast_name = {'dir color','D_para intra','D_para extra','D_orth extra','intra volfrac'};
          h = uicontextmenu;
          for k = 1:length(contrast_name)
             uimenu(h,'Label',contrast_name{k},'Callback',@(x,y) showFibs(ftr,k));
          end;
          set(gca,'uicontextmenu',h);
          
          ssds = get(fh,'userdata');
          ssds.fibcoloring = coloring;
          set(fh,'userdata',ssds);
          axis([0 25 0 25 0 3])
          view(0,90)

          try 
            aviobj = evalin('base','aviobj');
            F = getframe(gca);
            aviobj = addframe(aviobj,F);
            assignin('base','aviobj',aviobj);
            
          catch
              
          end;
          
          
          
          
function browsesig_Callback(h, eventdata, handles, varargin)      
          fh = findobj('tag','mesoGT_statfigure');
          ssds = get(fh,'userdata');              
          if not(isfield(ssds,'pos'));
              combobox_Callback([],[],[]) ; %% when called for first time 
              ssds = get(fh,'userdata');      
          end;          
          hr = mrstruct_read(ssds.ds.signal_fname);
          figure(999); stackview_compare(hr,ssds.ds.ftr)
          

function showMaps
          fh = findobj('tag','mesoGT_statfigure');
          ssds = get(fh,'userdata');              
          if not(isfield(ssds,'pos'));
              combobox_Callback([],[],[]) ; %% when called for first time 
              ssds = get(fh,'userdata');      
          end;          
          contrast_name = {'intra Dax','extra Dax','extra Drad','vf','vf_csf','Segment density','Terminal Density'};
          ran = {[1 3], [1 3], [0 2],[0 1],[0 0.2], [0 4], [0 4]};
          for j = 1:4,
              h{j} = uicontextmenu;
              for k = 1:length(contrast_name)
                 uimenu(h{j},'Label',contrast_name{k},'Callback',@(x,y) mapSel(x,y,j,k));
              end;
          end;          
          
          
          i = round(ssds.pos);
          
          if ssds.view == 3,
            x = flipdim(permute(squeeze(ssds.maps(:,:,i,:)),[2 1 3]),1);
            x = squeeze(ssds.maps(:,:,i,:));
          elseif ssds.view == 2,
            x = flipdim(permute(squeeze(ssds.maps(:,i,:,:)),[2 1 3]),1);
          elseif ssds.view == 1,
            x = flipdim(permute(squeeze(ssds.maps(i,:,:,:)),[2 1 3]),1);
          end;
          pos = [1 2 5 6];
          for k = 1:4,
              subplot(4,4,pos(k));
              hss = imagesc(x(:,:,ssds.contrast(k)),ran{ssds.contrast(k)}); title('intra axonal'); axis off;
%              hss = imagesc(x(:,:,ssds.contrast(k))); title('intra axonal'); axis off;
              title(contrast_name{ssds.contrast(k)});
              gcapos = get(gca,'position');
              set(gca,'position',[gcapos(1:2) gcapos(3:4)*1.25]);
              set(hss,'uicontextmenu',h{k});
          end;
          
          
function combobox_Callback(h, eventdata, handles, varargin)    
    combo = findobj(gcf,'tag','mesoGT_slicepopup');    
    fh = findobj(gcf,'tag','mesoGT_statfigure');
    ssds = get(fh,'userdata');
    ssds.view = get(combo,'value');
    
    slider = findobj('tag','mesoGT_slider');  
    slider = slider(1);
    set(slider,'min',1);
    set(slider,'max',size(ssds.maps,ssds.view(1)));
    set(slider,'value',round(size(ssds.maps,ssds.view(1))/2));
    ssds.pos = round(size(ssds.maps,ssds.view(1))/2);
    set(slider,'sliderstep',[0.01 0.1]);
    set(fh,'userdata',ssds);
    showMaps;

             
function slider_Callback(h, eventdata, handles, varargin)    
    slider = findobj(gcf,'tag','mesoGT_slider');    
    fh = findobj(gcf,'tag','mesoGT_statfigure');
    ssds = get(fh,'userdata');
    ssds.pos = get(slider,'value');
    set(fh,'userdata',ssds);
    showMaps;

    
function mapSel(x,y,j,k)
    fh = findobj(gcf,'tag','mesoGT_statfigure');
    ssds = get(fh,'userdata');
    ssds.contrast(j) = k;    
    set(fh,'userdata',ssds);
    showMaps;
              

function plotHists(datastruct,maps)
%%
        params = datastruct.params;
        vox = datastruct.vox;
        cD =  cat(2,datastruct.ftr.curveD{:});




        wn = 4; wm = 4;

    
       lw = 2;
       subplot(wn,wm,9);
       
       mask = not(isnan(maps(:,:,:,1)));
       mask = mask(:);
       
       b = 0:0.1:3;
       xx = maps(:,:,:,1); h1 = histc(xx(mask),b);
       xx = maps(:,:,:,2); h2 = histc(xx(mask),b);
       xx = maps(:,:,:,3); h3 = histc(xx(mask),b);
%        h1 = histc(cD(1,:),b);
%        h2 = histc(cD(2,:),b);
%        h3 = histc(cD(3,:),b);
       
       plot(b+0.01,h1/max(h1(:)),'g','linewidth',lw); hold on;
       plot(b+0.01,h2/max(h2(:)),'b','linewidth',lw); hold on;
       plot(b+0.01,h3/max(h3(:)),'b--','linewidth',lw); hold on;
       hold off;
       ax = axis; axis([0 3 ax(3) ax(4)]); grid on; set(gca,'yticklabel',[])       
       
       grid on; set(gca,'yticklabel',[])
       title('Diffusion parameters');
       
       subplot(wn,wm,10);
       b = 0:0.02:1;
       xx = maps(:,:,:,4); h1 = histc(xx(mask),b);
       %h1 = histc(cD(4,:),b);
       
       plot(b+0.01,h1/max(h1(:)),'g','linewidth',lw); hold on;
       hold off;
       ax = axis; axis([0 1 ax(3) ax(4)]); grid on; set(gca,'yticklabel',[])
       title('Volume fractions');
       
       
       subplot(wn,wm,11);
%        [h b] = myhist(cD(4,:),1+0*cD(4,:),[0 1],40);
%        plot(b*params.p_weight,h/max(h(:)),'r','linewidth',lw);  hold on;              
%        [h b] = myhist(cD(5,:),1+0*cD(4,:),[0 1],40);
%        plot(b*params.p_weight,h/max(h(:)),'b','linewidth',lw);  hold on;              
%        vfcorrwe = maps(:,:,:,4);
%        b2 = 0:0.05:1.5;
%        h2 = histc(vfcorrwe(vfcorrwe(:)>0),b2);       
%        plot(b2+0.025,h2/max(h2(:)),'g','linewidth',lw); 
%        hold off; grid on; set(gca,'yticklabel',[])
%         legend('w_D','w_q','S(0)');
%        title('weight')

%%
       subplot(wn,wm,12);   
       
       sz = size(datastruct.vfmap);
       osamp =  (datastruct.ftr.trackParam.params.p_wid); 
       pos =  datastruct.state(1:3,:)'*osamp;
       pos = pos ./ repmat(datastruct.ftr.vox(:)',[size(pos,1) 1]);
       mask = datastruct.vfmap(:,:,:,1)>0;
       
       fd = hist4(pos,sz(1:3));
       bar(0:6,histc(fd(mask(:)),-0.5:1:5.5));
       axis(axis.*[ 0 0 1 1] + [-1 7 0 0]);
       title('segments per voxel')
  %%     
       
%        fd = maps(:,:,:,5);
%        bar(0:6,histc(fd(maski(:)),-0.5:1:5.5));
%        axis(axis.*[ 0 0 1 1] + [-1 7 0 0]);
%        title('segments per voxel')

%        
%        b2 = 0:0.05:3;
%        b = 0:0.05:3;
%        subplot(wn,wm,13);
%        xh = hist3([cD(2,:)' cD(3,:)'],'edges',[{b2} {b}]);
%        imagesc(b,b2,xh(1:end-1,1:end-1));          
%        xlabel('Dextra perp'); ylabel('Dextra para');       
% 
%        subplot(wn,wm,14);
%        xh = hist3([cD(1,:)' cD(3,:)'],'edges',[{b2} {b}]);
%        imagesc(b,b2,xh(1:end-1,1:end-1));           
%        xlabel('Dextra perp'); ylabel('Dintra para');       
%        
%        subplot(wn,wm,15);
%        b = 0:0.015:1;
%        xh = hist3([cD(1,:)' cD(4,:)'],'edges',[{b2} {b}]);
%        imagesc(b,b2,xh(1:end-1,1:end-1)); hold on;    
%        xlabel('v_i'); ylabel('Dintra para');       
% 
%        subplot(wn,wm,16);
%        
%        cla
%        b2 = 0:0.015:1;
%        b = 0:0.015:1;
%        xh = hist3([cD(3,:)' cD(4,:)'],'edges',[{b2} {b}]);
%        imagesc(b,b2,xh(1:end-1,1:end-1)); hold on;
%        xlabel('v_i'); ylabel('Dextra perp');       
%%
function h = hist4(d,sz)
    d = double(floor(d)+1);

    d = d(d(:,1)>=1 & d(:,1)<=sz(1) & d(:,2)>=1 & d(:,2)<=sz(2) & d(:,3)>=1 & d(:,3)<=sz(3) ,:);

    didx = sub2ind(sz,d(:,1),d(:,2),d(:,3));
    h = sparse([didx(:);sz(1)*sz(2)*sz(3)],ones(length(didx(:))+1,1),ones(length(didx(:))+1,1));
    h = reshape(full(h),sz);




              