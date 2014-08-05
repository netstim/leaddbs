function BW=ls_bwareaopen(BW, P, CONN)
    CC=ls_bwconncomp(BW, CONN);
    for i=1:CC.NumObjects
        if (numel(CC.PixelIdxList{i}) < P)
            % disp(['Removing ' num2str(i)]);
            BW(CC.PixelIdxList{i}) = 0;
        end
    end
end