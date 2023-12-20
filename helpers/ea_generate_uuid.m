function uid = ea_generate_uuid
% Java built-in type 4 random UUID generator

uid = char(matlab.lang.internal.uuid());
