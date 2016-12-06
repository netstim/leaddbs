function ea_checkspm


ver=spm('version');


if str2double(ver(8:11))<6906
    msgbox('Some functions (such as SPM SHOOT and DARTEL) may not be available using your SPM version. Please upgrade SPM12 to at least revision 6906 (or simply update to newest release).');
end