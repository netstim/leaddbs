function id=ea_method2methodid(obj)
switch obj.statmetric
    case 1
        id='ttests';
    case 2
        id='spearman';
        switch obj.efieldmetric
            case 'mean'
                id=[id,'_mean'];
            case 'peak'
                id=[id,'_peak'];
            case 'sum'
                id=[id,'_sum'];
        end
end
end