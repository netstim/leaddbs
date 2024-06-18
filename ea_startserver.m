function ea_startserver(handles, bids)
    [jsonData] = ea_getTutorData(handles, bids); 
    [mac_address] = ea_getMac();
%     data = struct('key1', 'value1', 'key2', 123, 'key3', [1, 2, 3]);
%     jsonData = jsonencode(data); % Will need to replace this with the actual user data
%     baseUrl = 'http://localhost:3000/?userid=';
    baseUrl = 'http://87.106.232.66/?userid=';
    url = 'http://87.106.232.66:8000/my-endpoint';  % URL of your FastAPI endpoint
    options = weboptions('RequestMethod', 'post', 'MediaType', 'application/json', 'HeaderFields', {'User-ID', mac_address});
    response = webwrite(url, jsonData, options);
%     disp(response)
    userUrl = strcat(baseUrl, mac_address);
    web(userUrl, '-browser');
end
% data = struct('key1', 'value1', 'key2', 123, 'key3', [1, 2, 3]);
% jsonData = jsonencode(data);
% url = 'http://localhost:8000/my-endpoint';  % URL of your FastAPI endpoint
% options = weboptions('RequestMethod', 'post', 'MediaType', 'application/json');
% response = webwrite(url, jsonData, options);
% url = 'http://localhost:8000/my-endpoint';
% options = weboptions('RequestMethod', 'post', 'MediaType', 'application/json');
% response = webwrite(url, jsonData, options);
% url = 'http://localhost:8000/'
% response = webwrite(url, jsonData, options);
% url = 'http://localhost:8000/';  % URL of your FastAPI endpoint
% data = webread(url);
% disp(data);
% url = 'http://localhost:8000/my-endpoint';  % URL of your FastAPI endpoint
% response = webwrite(url, jsonData, options);
% [status, fullname] = system('id -F')
% [s, out] = dos('vol');
% networkinterfaces = java.net.NetworkInterface.getNetworkInterfaces
% nimacs = cell(0, 2);
% while networkinterfaces.hasMoreElements
% networkinterface = networkinterfaces.nextElement;
% macstring = strjoin(cellstr(dec2hex(typecast(networkinterface.getHardwareAddress, 'uint8'))), ':');
% nimacs = [nimacs; {char(networkinterface.getDisplayName), macstring}];
% end
% cell2table(nimacs, 'VariableNames', {'Interface', 'MAC'})
% vol
% ea_getMac
% url2 = 'http://localhost:8000/items/'
% options = weboptions('HeaderFields', {'User-Agent', headers.User-Agent});
% options = weboptions('HeaderFields', {'User-ID', mac_adress});
% options = weboptions('HeaderFields', {'User-ID', mac_address});
% response = webread(url, options);
% response = webread(url2, options);
% response
% response = webread(url2, options);
% response
% options = weboptions('HeaderFields', {'User-ID', mac_address});
% response = webread(url2, options);
% response
% response = webread(url2, options);
% response = webwrite(url, jsonData, options);
% mac_address