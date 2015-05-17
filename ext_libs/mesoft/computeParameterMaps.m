
function maps = computeParameterMaps(ftr,M);

      

  
      idxmask = ftr.user.S2(:,:,:,2)>0;
      
      % kernel for gathering values withtin original DWI-voxel
      osamp =  (ftr.trackParam.params.p_wid); 
      osamp = osamp/M;
      st = max(ceil(osamp/2),1);
      ker = ones(osamp,osamp,osamp); %ker = ker /sum(ker(:));
      
      % compute averages
      %  vf(:,:,:,1) - raw vfi 
      %  vf(:,:,:,2) - projection of raw model onto constant
      %  vf(:,:,:,3) - segment's contributions within voxel  \sum_i w_i
      vf = ftr.user.vf; vf(vf==1)=0;
      vfsm = imfilter(vf,ker); % vfsm(vfsm==1) = 0;
      
      % get lowres map with #voxels in highres map
      acpervox = imfilter(double(vf(:,:,:,1)>0),ker); 
      acpervox = acpervox(st:osamp:end,st:osamp:end,st:osamp:end,:);
      
      % downsample vfi to original resolution
      vfsm = vfsm(st:osamp:end,st:osamp:end,st:osamp:end,:);
      vfsm = vfsm ./ (eps+repmat(acpervox,[1 1 1 3]));
      
      % mean of signal 
      vfmsig = zeros(size(ftr.user.b0avg)); 
      vfmsig(idxmask) = ftr.user.meansignal;
      sz = size(vfsm);
      if any(size(vfmsig) ~= sz(1:3)),
          fac = sz(1)/size(vfmsig,1);
          tmp = vfmsig;
          vfmsig = zeros(sz(1:3));
          for k = 1:fac,
              for j = 1:fac,
                  for i = 1:fac,
                    vfmsig(k:fac:end,j:fac:end,i:fac:end) = tmp;
                  end;
              end;
          end;          
      end;
      
      if  ftr.trackParam.params.alpha == 1, %% with implicit SW

          % compute remaining meansignal without vfsw
          vfsm(:,:,:,2) = vfmsig-vfsm(:,:,:,2);
          vfsm2tmp = vfsm(:,:,:,2); 
          vfsm2tmp(vfsm2tmp<0) = 0; 
          vfsm(:,:,:,2) = vfsm2tmp;

          % compute predicted S0
          S0 = (vfsm(:,:,:,3)+vfsm(:,:,:,2));

          % compute vfi normalized by S0
          vfsm(:,:,:,1) = (vfsm(:,:,:,1) .* vfsm(:,:,:,3)) ./ S0;

          % compute vfsw normalized by S0
          vfsm(:,:,:,2) = vfsm(:,:,:,2) ./ S0;

          % compute vfe
          vfsm(:,:,:,3) = 1-sum(vfsm(:,:,:,1:2),4);
          
          % thats' S0
          vfsm(:,:,:,4) = S0;
      else
          S0 = vfsm(:,:,:,3);
          vfsm(:,:,:,1) = vfsm(:,:,:,1);
          vfsm(:,:,:,2) = vfsm(:,:,:,1)*0;
          vfsm(:,:,:,3) = (1-vfsm(:,:,:,1)).*(vfsm(:,:,:,1)>0);
          vfsm(:,:,:,4) = S0;
      end;
          
      
      % mask out with nans
      maski = vfsm(:,:,:,1)>0;
      for k = 1:4,
          tmp = vfsm(:,:,:,k);
          tmp(not(maski)) = nan;
          vfsm(:,:,:,k) = tmp;
      end;
      

      % compute average segment count in voxel
      sz = size(vfsm);
      pos = ftr.user.P(1:3,:)';
      pos = pos ./ repmat(ftr.vox(:)',[size(pos,1) 1]);
      fdsm = hist4(pos,sz(1:3));
      fdsm = fdsm ./ (acpervox+eps);
            
      
      % compute average endpoint count in voxel
      idx = ftr.user.P(9,:) == -1 | ftr.user.P(10,:) == -1;
      ep = hist4(pos(idx,:),sz(1:3));      
      ep = ep ./ (acpervox+eps);

      
 
      maps.vfi = vfsm(:,:,:,1);
      maps.vfsw = vfsm(:,:,:,2);
      maps.vfe = vfsm(:,:,:,3);
      maps.S0 = vfsm(:,:,:,4);
      maps.segcount = fdsm;
      maps.termcount = ep;
      
      
      
      
      
function h = hist4(d,sz)
d = double(floor(d)+1);

d = d(d(:,1)>=1 & d(:,1)<=sz(1) & d(:,2)>=1 & d(:,2)<=sz(2) & d(:,3)>=1 & d(:,3)<=sz(3) ,:);

didx = sub2ind(sz,d(:,1),d(:,2),d(:,3));
h = sparse([didx(:);sz(1)*sz(2)*sz(3)],ones(length(didx(:))+1,1),ones(length(didx(:))+1,1));
h = reshape(full(h),sz);

      


      
      
      
      
      