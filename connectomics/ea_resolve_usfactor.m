function factor=ea_resolve_usfactor(upsample)

switch upsample.factor
    case 1
        factor=1;
    case 2
        factor=1.5;
    case 3
        factor=2;
    case 4
        factor=3;
    case 5
        factor=4;
end