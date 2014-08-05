function CC=ls_bwconncomp(BW, CONN)
    switch CONN
        case 4
            CONN=6;
        case 8
            CONN=18;
    end
        
    [L, NUM] = spm_bwlabel(double(BW), CONN);
    
    CC.Connectivity = CONN;
    CC.ImageSize = size(BW);
    CC.NumObjects = NUM;
    CC.PixelIdxList = ls_conncomp_pix_list(L, NUM);
end
