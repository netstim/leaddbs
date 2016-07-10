function uid = generate_guid
%GENERATE_GUID  Generate a globally unique ID.
%   UID = GENERATE_GUID creates a universally unique identifier (UUID). 
%   A UUID represents a 128-bit value. For more information including 
%   algorithms used to create UUIDs, see RFC 4122: 
%   A Universally Unique IDentifier (UUID) URN Namespace, 
%   section 4.2 "Algorithms for Creating a Time-Based UUID".

%   Copyright 2012 Changjiang Yang
%   my_last_name.cj@gmail.com

import java.util.UUID;

uid = char(UUID.randomUUID());

return