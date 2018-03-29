function pred=ea_predict_Bmap(bmap,connfile,sk,usemask)
if ~exist('sk','var')
    sk='';
end
if ~exist('usemask','var');
    usemask='';
end


[pth,fn,ext]=fileparts(bmap);
b1=ea_load_nii(bmap);
intercept=ea_load_nii(fullfile(pth,[fn,'_intercept',ext]));
if exist(fullfile(pth,[fn,'_deviation',ext]),'file')
    dev=ea_load_nii(fullfile(pth,[fn,'_deviation',ext]));
    allzeros=(b1.img==0); % add all zeros to zeromask
    allzeros=logical(allzeros+isnan(b1.img)); % add all nans to zeromask
    
    weights=1./(8*dev.img);
    weights=weights-ea_nanmin(weights(:)); % get rid of negatives
    weights(allzeros)=0; % set everything to zero that was all zeros or nans
    weights=weights./sum(abs(weights(:)));
else % robustfit
    weights=b1.img(:)./sum(b1.img(:));
end


ptconn=ea_genX(connfile,1,'tmp.nii',ea_getmask(usemask),'');

pred=intercept.img(:)+b1.img(:).*ptconn(:);

pred=ea_nansum(pred.*weights(:));
    
    