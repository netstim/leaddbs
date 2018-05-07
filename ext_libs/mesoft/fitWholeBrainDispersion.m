function results = fitWholeBrainDispersion(hr)
    


    ds.maps =  compMesoParams(hr);
    ds.hr = hr;
    dirs = load('dwidirections');
    ds.dirs256 = [dirs.dirs256 -dirs.dirs256];
    clear dirs;

    modelstruc = evalin('base','modelstruc');
    
    proto.ten = hr.user.bTensor;
    proto.b = squeeze(round((proto.ten(1,1,:)+proto.ten(2,2,:)+proto.ten(3,3,:))));     
    proto.buni = unique(round(proto.b/100))*100;  
    for k = 1:size(proto.ten,3);
        [U D] = eigs(proto.ten(:,:,k));
        proto.dirs(:,k) = U(:,1);
    end;
    
    
    
        
    
    
    for x = 1:size(hr.dataAy,1),
    for y = 1:size(hr.dataAy,2),
    for z = 1:size(hr.dataAy,3),

        fprintf('%i %i %i / %i %i %i \n',x,y,z,size(hr.dataAy,1),size(hr.dataAy,2),size(hr.dataAy,3));

        c = [x y z];

        Smeas = squeeze(ds.hr.dataAy(c(1),c(2),c(3),:));
        Smeas = Smeas ./ mean(Smeas(round(proto.b/100)==0));

        
        %% get parameters of segment
                        
        Dpara_in =  ds.maps(c(1),c(2),c(3),1);
        Dpara_ex =  ds.maps(c(1),c(2),c(3),2)+ds.maps(c(1),c(2),c(3),3);
        Dorth =  ds.maps(c(1),c(2),c(3),3);
        vf =        ds.maps(c(1),c(2),c(3),4);
        vf(vf>1) = 1;
        vf_csf = ds.maps(c(1),c(2),c(3),5);
        vf_csf(vf_csf<0) = 0;
        snr = ds.maps(c(1),c(2),c(3),7);
        
        
        
        
        disptypes = {'watson','poisson','mises','heat'};
        res = fitDispersion(Dpara_in,Dpara_ex,Dorth,vf,vf_csf,snr,modelstruc.noisedeg,proto, Smeas,ds.dirs256,disptypes)
        res = rmfield(res,'S');
        res = rmfield(res,'fod');
        
        if not(exist('results')),
            results = repmat(res,[size(hr.dataAy,1) size(hr.dataAy,2) size(hr.dataAy,3)]);
        end;
        results(x,y,z) = res;
        
        
        
    end
    end
    end;
