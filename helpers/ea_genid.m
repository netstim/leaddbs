function ID = ea_genid(ID)
% Generate patient ID based on MD5 hashtag.

Hash = GetMD5(ID);
ID = Hash(1:5);
