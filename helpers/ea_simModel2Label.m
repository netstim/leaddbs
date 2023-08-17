function label = ea_simModel2Label(simModel)
% Return a label for specific simulation model (of VTA).

switch simModel
    case 'Dembek 2017'
        label = 'dembek';
    case 'Fastfield (Baniasadi 2020)'
        label = 'fastfield';
    case 'Kuncel 2008'
        label = 'kuncel';
    case 'Maedler 2012'
        label = 'maedler';
    case 'OSS-DBS (Butenko 2020)'
        label = 'ossdbs';
    case 'SimBio/FieldTrip (see Horn 2017)'
        label = 'simbio';
    otherwise
        error('Simulation method not recognizable!')
end
