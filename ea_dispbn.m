function ea_dispbn(varargin)
if nargin>0
    if strcmp(varargin{1},'ee');
        
        for space=1:100
            clc
            spstr=repmat(' ',1,space);
            %% display banner..
            
            disp('                 Welcome to                            ');
            disp('                                                            ');
            disp('  _______________________________________                   ');
            disp(' /           ______     ___      _____   \                  ');
            disp([' | | ocali-  |         /   \    |     \  | ',spstr,'       _______   ']);
            disp([' | | zation  | lec     |   |    | BS   \ | ',spstr,'    = / O   / \  ']);
            disp([' | |         |---      /___\    | Sur- | | ',spstr,'   = (   \./___) ']);
            disp([' | | of      | trodes |     |   | gery / | ',spstr,'  =   \_ __ _)   ']);
            disp([' | |_____    |_____   / fter \  |_____/  | ',spstr,'   =   (O)(O)    ']);
            disp(' \_______________________________________/                  ');
            disp('                                                            ');
            disp('                   Toolbox.                            ');
            disp(['                     ',ea_getvsn,'                    ']);
            pause(0.001)
        end
        
        
        clc
        
        disp('                 Welcome to                ');
        disp('                                           ');
        disp('  _______________________________________  ');
        disp(' /           ______     ___      _____   \ ');
        disp(' | | ocali-  |         /   \    |     \  | ');
        disp(' | | zation  | lec     |   |    | BS   \ | ');
        disp(' | |         |---      /___\    | Sur- | | ');
        disp(' | | of      | trodes |     |   | gery / | ');
        disp(' | |_____    |_____   / fter \  |_____/  | ');
        disp(' \_______________________________________/ ');
        disp('                                           ');
        disp('                   Toolbox.                ');
        disp(['                     ',ea_getvsn,'        ']);
    end
else
    space=5;
    clc
    spstr=repmat(' ',1,space);
    %% display banner..
    
    disp('                 Welcome to                            ');
    disp('                                                            ');
    disp('  _______________________________________                   ');
    disp(' /           ______     ___      _____   \                  ');
    disp([' | | ocali-  |         /   \    |     \  | ',spstr,'       _______   ']);
    disp([' | | zation  | lec     |   |    | BS   \ | ',spstr,'    = / O   / \  ']);
    disp([' | |         |---      /___\    | Sur- | | ',spstr,'   = (   \./___) ']);
    disp([' | | of      | trodes |     |   | gery / | ',spstr,'  =   \_ __ _)   ']);
    disp([' | |_____    |_____   / fter \  |_____/  | ',spstr,'   =   (O)(O)    ']);
    disp(' \_______________________________________/                  ');
    disp('                                                            ');
    disp('                   Toolbox.                            ');
    disp(['                     ',ea_getvsn,'                    ']);
    pause(0.001)
    
    
end