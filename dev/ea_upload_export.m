function response=ea_upload_export(directory,patientPseudonym,usercredentials)

import matlab.net.http.io.FileProvider
import matlab.net.http.RequestMessage
import matlab.net.http.HeaderField
import matlab.net.http.io.MultipartProvider
import matlab.net.http.io.MultipartFormProvider
import matlab.net.http.io.FileProvider
import matlab.net.http.MessageBody
import matlab.net.http.field.AuthorizationField
import matlab.net.http.Credentials
import matlab.net.http.HTTPOptions
import matlab.net.http.MediaType

options=ea_getptopts(directory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input from user and LEAD-DBS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filePath = [options.root,options.patientname,filesep,'export',filesep,'zip',filesep,options.patientname,'.zip'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Disable internet proxy in Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code to automatically refresh token (in progress ...)
%{ 
input = struct('username',options.prefs.tbase.user,'password',options.prefs.tbase.pw, 'grant_type', 'password', 'client_id', options.prefs.tbase.clientId, 'client_secret', options.prefs.tbase.clientSecret, 'response_type', 'token');
inputParameters = jsonencode(input);
body = MessageBody(inputParameters);
headerField = matlab.net.http.field.ContentTypeField('application/x-www-form-urlencoded');
tokenReq = RequestMessage('post',headerField,body);
show(tokenReq);
response = tokenReq.send(options.prefs.tbase.authUrl)
%}
prefs=ea_prefs;

user = usercredentials{1}; % to be replaced by technical user
pw = usercredentials{2};


% Obtain token
% response = webwrite(prefs.tbase.authUrl,'username',user,'password',pw,'client_id',prefs.tbase.clientId,'client_secret',prefs.tbase.clientSecret,'grant_type','password');
% token = ['Bearer ' response.access_token];
o = weboptions('CertificateFilename',prefs.tbase.certificateFile);
response = webwrite(prefs.tbase.authUrl,'username',user,'password',pw,'client_id',prefs.tbase.clientId,'client_secret',prefs.tbase.clientSecret,'grant_type','password',o);
token = ['Bearer ' response.access_token];

% Send request
uploadUrl = prefs.tbase.endpoint + patientPseudonym;
fp = FileProvider(filePath);

formProvider = MultipartFormProvider("file",fp);
headerFields = HeaderField('Authorization', token,'Content-Type','multipart/form-data');
options = matlab.net.http.HTTPOptions('ConnectTimeout',20,'CertificateFilename',prefs.tbase.certificateFile);
req = RequestMessage('post',headerFields,formProvider);
response = req.send(uploadUrl,options);

% Send request
%{ 

uploadUrl = prefs.tbase.endpoint + patientPseudonym;
fp = FileProvider(filePath);
formProvider = MultipartFormProvider("file",fp);
headerFields = HeaderField('Authorization', token,'Content-Type','multipart/form-data');
options = matlab.net.http.HTTPOptions('ConnectTimeout',20);
req = RequestMessage('post',headerFields,formProvider);
response = req.send(uploadUrl,options);
%}
