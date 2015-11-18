function newmarkers=ea_correctcoords(oldmarkers,trajectory,command)
newmarkers=oldmarkers;

grone=[0.5,3]; % step-sizes

switch command.Key
    case 'downarrow'
        % side=1
        
        
        try
            
            trajvector=gettraj(trajectory,1);
            side=1;
            newmarkers(side).head=oldmarkers(side).head+grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).tail=oldmarkers(side).tail+grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).x=oldmarkers(side).x+grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).y=oldmarkers(side).y+grone(1+ismember('shift',command.Modifier))*trajvector;
            
        end

    case 'uparrow'
        % side=1
        try
            trajvector=gettraj(trajectory,1);
            side=1;
            newmarkers(side).head=oldmarkers(side).head-grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).tail=oldmarkers(side).tail-grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).x=oldmarkers(side).x-grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).y=oldmarkers(side).y-grone(1+ismember('shift',command.Modifier))*trajvector;

        end

    case 'rightarrow'
        % side=2
        try
            
            trajvector=gettraj(trajectory,2);
            side=2;
            newmarkers(side).head=oldmarkers(side).head+grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).tail=oldmarkers(side).tail+grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).x=oldmarkers(side).x+grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).y=oldmarkers(side).y+grone(1+ismember('shift',command.Modifier))*trajvector;

        end

    case 'leftarrow'
        % side=2
        try
            trajvector=gettraj(trajectory,2);
            side=2;
            newmarkers(side).head=oldmarkers(side).head-grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).tail=oldmarkers(side).tail-grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).x=oldmarkers(side).x-grone(1+ismember('shift',command.Modifier))*trajvector;
            newmarkers(side).y=oldmarkers(side).y-grone(1+ismember('shift',command.Modifier))*trajvector;

        end

        
end
switch command.Character
    case {'-','_'}
        try
            for side=1:2
               trajvector=gettraj(trajectory,side);
                newmarkers(side).tail=oldmarkers(side).tail+trajvector.*(0.3*(ismember('shift',command.Modifier)+1));
            end
            
        end
        
    case {'+','*'}
        try
            
            for side=1:2
              trajvector=gettraj(trajectory,side);
              newmarkers(side).tail=oldmarkers(side).tail-trajvector.*(0.3*(ismember('shift',command.Modifier)+1));
            end
        end
        
    case 's'
        %newmarkers(side).x=rotate(
        
        
        
end
%disp(command);



function traj=gettraj(trajectory,side)

traj=mean(diff(trajectory{side}));