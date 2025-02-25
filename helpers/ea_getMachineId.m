function machineId = ea_getMachineId
% Get unique machine ID

if ispc
    machineId = getPCUUID;

    if isempty(machineId)
        machineId = convertToUUID(getPCCPUId);
        if isempty(machineId)
            machineId = convertToUUID(getPCBIOSSerial);
        end
    end
elseif ismac
    [~, cmdout] = system('system_profiler SPHardwareDataType | grep "Hardware UUID" | awk ''{print $3}''');
    machineId = '';%strip(cmdout);

    if isempty(machineId)
        [~, cmdout] = system('system_profiler SPHardwareDataType | grep "Serial Number (system)" | awk ''{print $4}''');
        machineId = convertToUUID(strip(cmdout));
    end
elseif unix
    [~, cmdout] = system('cat /etc/machine-id 2>/dev/null');
    machineId = convertToUUID(strip(cmdout));

    if isempty(machineId)
        [~, cmdout] = system('cat /var/lib/dbus/machine-id 2>/dev/null');
        machineId = convertToUUID(strip(cmdout));
    end
end


function uuid = getPCUUID
[~, cmdout] = system('wmic csproduct get UUID');
cmdout = strsplit(strip(cmdout), '\n');
if numel(cmdout) >= 2
    uuid = cmdout{2};
else
    uuid = '';
end


function cpu_id = getPCCPUId
[~, cmdout] = system('wmic cpu get ProcessorId');
cmdout = strsplit(strip(cmdout), '\n');
if numel(cmdout) >= 2
    cpu_id = cmdout{2};
else
    cpu_id = '';
end


function bios_id = getPCBIOSSerial
[~, cmdout] = system('wmic bios get serialnumber');
cmdout = strsplit(strip(cmdout), '\n');
if numel(cmdout) >= 2
    bios_id = cmdout{2};
else
    bios_id = '';
end


function uuid = convertToUUID(str)
if isempty(str)
    uuid = '';
    return;
end

if isUUIDFormat(str)
    uuid = str;
    return;
end

hash = str2hash(str);

uuid = sprintf('%08X-%04X-%04X-%04X-%012X', ...
    hash(1), ...
    bitand(hash(2), hex2dec('FFFF')), ...
    bitand(bitor(hex2dec('4000'), bitand(hash(3), hex2dec('0FFF'))), hex2dec('FFFF')), ... % Version 4 UUID
    bitand(bitor(hex2dec('8000'), bitand(hash(4), hex2dec('3FFF'))), hex2dec('FFFF')), ... % Variant 1 UUID
    hash(5));


function hash = str2hash(str)
import java.security.*;
import java.math.*;

md = MessageDigest.getInstance('SHA-256');
hashBytes = md.digest(uint8(str));
hash = typecast(hashBytes(1:20), 'uint32');


function isUUID = isUUIDFormat(str)
pattern = '^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$';
isUUID = ~isempty(regexp(str, pattern, 'once'));
