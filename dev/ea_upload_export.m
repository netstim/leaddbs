function response=ea_upload_export(directory)

options=ea_getptopts(directory);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input from user and LEAD-DBS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([directory,filesep,'ea_pseudonym.mat'],'file')
    patientPseudonym = inputdlg('Please enter patient synonym','Patient Synonym',1,{options.patientname});
    save([directory,filesep,'ea_pseudonym.mat'],'patientPseudonym');
else % if patient has already been exported, we already know the pseudonym.
    p=load([directory,filesep,'ea_pseudonym.mat']);
    patientPseudonym=p.patientPseudonym; clear p
end
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

answ = inputdlg({'Username','Password'},'Enter User Credentials',[1 35],{prefs.tbase.user,prefs.tbase.pw});
user = answ{1}; % to be replaced by technical user
pw = answ{2};


% Obtain token
response = webwrite(prefs.tbase.authUrl,'username',user,'password',pw,'client_id',prefs.tbase.clientId,'client_secret',prefs.tbase.clientSecret,'grant_type','password');
token = ['Bearer ' response.access_token];

% Send request
uploadUrl = prefs.tbase.endpoint + patientPseudonym;
fp = FileProvider(filePath);
formProvider = MultipartFormProvider("file",fp);
headerFields = HeaderField('Authorization', token,'Content-Type','multipart/form-data');
options = matlab.net.http.HTTPOptions('ConnectTimeout',20);
req = RequestMessage('post',headerFields,formProvider);
response = req.send(uploadUrl,options);
