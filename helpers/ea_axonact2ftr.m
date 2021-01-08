function ea_axonact2ftr(axonacts)
% Convert axon activation from OSS-DBS to FTR (and additional TRK) format.
%
% Input can be the folder where the activation CSV files are stored or the
% CSV file itself.


if isfolder(axonacts)
    axonacts = ea_regexpdir(axonacts, '^Activation', 0);
elseif isfile(axonacts)
    axonacts = {axonacts};
end

ftr.ea_fibformat = '1.0';
ftr.fourindex = 1;
ftr.voxmm = 'mm';

for f=1:length(axonacts)
    axonLength = regexp(axonacts{f}, '.*_(\d+)\.csv$', 'tokens', 'once');
    axonLength = str2double(axonLength{1});

    ftr.fibers = cell2mat(textscan(fopen(axonacts{f}), '%f %f %f %f', 'CollectOutput', 1));
    axonNum = size(ftr.fibers,1)/axonLength;

    ftr.fibers(:,4) = repelem(1:axonNum, axonLength)';
    ftr.idx = ones(axonNum,1)*axonLength;
    save(strrep(axonacts{f}, '.csv', '.mat'), '-struct', 'ftr');
    ea_ftr2trk(strrep(axonacts{f}, '.csv', '.mat'));
end
