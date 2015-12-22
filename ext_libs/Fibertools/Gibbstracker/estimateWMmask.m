function estimate = estimateWMmask(signal,tensor)

thres1 = 9/27;
thres2 = 0.25;
t3 = 0.8;
% thres1 = 9/27;
% thres2 = 0.25;
% t3 = 1;
% 
% thres1 = 16/27;
% thres2 = 0.5;
        idxb0 = find(squeeze(squeeze(tensor(1,1,:)+tensor(2,2,:)+tensor(3,3,:))) < 101);
        idxb1 = setdiff(1:size(tensor,3),idxb0);
        meanDWI = squeeze(mean(signal(:,:,:,idxb1),4));
        meanB0 = squeeze(mean(signal(:,:,:,idxb0),4));
        maxmeanDWI = max(meanDWI(:));   
        binsz = 0.005;
        h = histc(meanDWI(:)/maxmeanDWI,0:binsz:1);
        hs = real(ifft(fft(log(1+h)).*fftshift(fspecial('gaussian',[size(h,1) 1],30))));
        idx = find((hs > imdilate(hs,[1 0 1]')));
        [mm midx] = sort(hs(idx),'descend');
        thres = t3*(idx(midx(1)) + idx(midx(2)))/2*binsz*maxmeanDWI;
        estimate = ((meanDWI./meanB0 > thres2) .* (meanDWI>thres));
        estimate = smooth3(estimate,'box',3)>thres1;