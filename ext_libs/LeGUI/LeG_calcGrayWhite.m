function [ElecType,ElecIdxLin] = LeG_calcGrayWhite(app,ElecFullIdx)

%Calculates gray/white class for an electrode given the index location
%"ElecFullIdx". ElecFullIdx is a matrix of index values specifying the
%volume of the electrode.
%
%Tyler Davis
%20190629

Thresh = app.graywhite.threshold; %0.1
DefaultType = app.graywhite.defaulttype; %'Gray'

ElecIdxLin = sub2ind(size(app.GrayImg),ElecFullIdx(:,1),ElecFullIdx(:,2),ElecFullIdx(:,3));
GrayVals = app.GrayImg(ElecIdxLin);
WhiteVals = app.WhiteImg(ElecIdxLin);
mGrayVal = mean(GrayVals);
mWhiteVal = mean(WhiteVals);
if mGrayVal>Thresh || mWhiteVal>Thresh
    [~,h,stats] = ranksum(GrayVals,WhiteVals,'alpha',0.001);
    if h
        if isfield(stats,'zval')
            if stats.zval>0
                ElecType = 'Gray';
            else
                ElecType = 'White';
            end
        else
            if mGrayVal>mWhiteVal
                ElecType = 'Gray';
            else
                ElecType = 'White';
            end
        end
    else %distributions are not different but at least one type is >0.1
        switch stats.ranksum
            case 1 %if only one voxel provided (i.e. click on 2D away from electrode)
                ElecType = 'White';
            case 2
                ElecType = 'Gray';
            otherwise
                ElecType = DefaultType; %'Gray'
        end
    end
else
    ElecType = 'Unknown';
end