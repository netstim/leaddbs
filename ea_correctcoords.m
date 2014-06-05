function newcoords=ea_correctcoords(oldcoords,trajectory,command)
newcoords=oldcoords;

grone=[0.5,3]; % step-sizes

switch command.Key
    case 'downarrow'
        % side=1
        
        
        try
            
            trajvector=gettraj(trajectory,1);
            side=1;
            for el=1:4
                newcoords{side}(el,:)=oldcoords{side}(el,:)+grone(1+ismember('shift',command.Modifier))*trajvector;
            end
        end

    case 'uparrow'
        % side=1
        try
            trajvector=gettraj(trajectory,1);
            side=1;
            for el=1:4
                newcoords{side}(el,:)=oldcoords{side}(el,:)-grone(1+ismember('shift',command.Modifier))*trajvector;
            end
        end

    case 'rightarrow'
        % side=2
        try
            
            trajvector=gettraj(trajectory,2);
            side=2;
            for el=1:4
                newcoords{side}(el,:)=oldcoords{side}(el,:)+grone(1+ismember('shift',command.Modifier))*trajvector;
            end
        end

    case 'leftarrow'
        % side=2
        try
            trajvector=gettraj(trajectory,2);
            side=2;
            for el=1:4
                newcoords{side}(el,:)=oldcoords{side}(el,:)-grone(1+ismember('shift',command.Modifier))*trajvector;
            end
        end

        
end
switch command.Character
    case {'-','_'}
        try
            trajvector=gettraj(trajectory,2);
            feather=[0,0.1,0.2,0.3];
            for side=1:2
                for el=1:4
                    newcoords{side}(el,:)=oldcoords{side}(el,:)+trajvector.*(feather(el)*(ismember('shift',command.Modifier)+1));
                end
            end
            
        end
        
    case {'+','*'}
        try
            feather=[0,0.1,0.2,0.3];
            
            trajvector=gettraj(trajectory,2);
            for side=1:2
                for el=1:4
                    newcoords{side}(el,:)=oldcoords{side}(el,:)-trajvector.*(feather(el)*(ismember('shift',command.Modifier)+1));
                end
            end
        end
        
end
%disp(command);



function traj=gettraj(trajectory,side)

traj=mean(diff(trajectory{side}));