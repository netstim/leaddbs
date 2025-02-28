function ID = ea_genid_rand(num, width, mode)
% Generate unique random patient ID
arguments
    num {mustBeNumeric, mustBeGreaterThan(num, 0)} = 1 % number of ID to be generated
    width {mustBeNumeric} = 4 % width of ID
    mode {mustBeMember(mode, {'char', 'number', 'mixed'})} = 'char'
end

switch mode
    case 'char'
        range = 65:90;
    case 'number'
        range = 48:57;
    case 'mixed'
        range = [48:57,65:90];  
end

if length(range)^width < num
    ea_cprintf('CmdWinErrors', 'Possible combinations less than desired number of ID!\n');
    return;
end

genid = @() char(range(randi([1 length(range)], 1, width)));

if num == 1
    ID = genid();
else
    ID = arrayfun(@(x) genid(), 1:num, 'Uni', 0)';
    while length(unique(ID)) ~= num
        ID = arrayfun(@(x) genid(), 1:num, 'Uni', 0)';
    end
end
