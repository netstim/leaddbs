function [planedim,onedim, secdim, dstring, lstring, Ltxt, Rtxt,plusminusc,plusminusr,plusminusl]=ea_getdims(tracor,side)

switch tracor
    
    case 1 % transversal images
        onedim=1;
        secdim=2;
        planedim=3;
        dstring='tra';
        lstring='z = ';
        Ltxt='M';
        Rtxt='R';
        plusminusc='plus';
        switch side
            case 1
                plusminusr='minus';
                plusminusl='plus';
            case 2
                plusminusr='plus';
                plusminusl='minus';
        end
    case 2 % coronal images
        onedim=1;
        secdim=3;
        planedim=2;
        dstring='cor';
        lstring='y = ';
        Ltxt='M';
        Rtxt='R';
        plusminusc='minus';
        
        switch side
            case 1
                plusminusr='minus';
                plusminusl='plus';
            case 2
                plusminusr='plus';
                plusminusl='minus';
        end
    case 3 % saggital images
        onedim=2;
        secdim=3;
        planedim=1;
        dstring='sag';
        lstring='x = ';
        Ltxt='P';
        Rtxt='A';
        plusminusc='minus';
        switch side
            case 1
                plusminusr='plus';
                plusminusl='minus';
            case 2
                plusminusr='plus';
                plusminusl='minus';
        end
end
