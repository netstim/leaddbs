function newmarkers=ea_correctcoords(oldmarkers,trajectory,command,options)
newmarkers=oldmarkers;

grone=[0.1,0.5,3]; % step-sizes

switch command.Key
    case 'uparrow'
        try
            trajvector=gettraj(trajectory,options.elside);
            newmarkers(options.elside).head=oldmarkers(options.elside).head+grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier))*trajvector;
            newmarkers(options.elside).tail=oldmarkers(options.elside).tail+grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier))*trajvector;
            newmarkers(options.elside).x=oldmarkers(options.elside).x+grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier))*trajvector;
            newmarkers(options.elside).y=oldmarkers(options.elside).y+grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier))*trajvector;
        end

    case 'downarrow'
        try
            trajvector=gettraj(trajectory,options.elside);
            newmarkers(options.elside).head=oldmarkers(options.elside).head-grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier))*trajvector;
            newmarkers(options.elside).tail=oldmarkers(options.elside).tail-grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier))*trajvector;
            newmarkers(options.elside).x=oldmarkers(options.elside).x-grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier))*trajvector;
            newmarkers(options.elside).y=oldmarkers(options.elside).y-grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier))*trajvector;
        end

%     case 'leftarrow'
%         % side=2
%         try
%
%             trajvector=gettraj(trajectory,2);
%             side=2;
%             newmarkers(side).head=oldmarkers(side).head+grone(1+ismember('shift',command.Modifier))*trajvector;
%             newmarkers(side).tail=oldmarkers(side).tail+grone(1+ismember('shift',command.Modifier))*trajvector;
%             newmarkers(side).x=oldmarkers(side).x+grone(1+ismember('shift',command.Modifier))*trajvector;
%             newmarkers(side).y=oldmarkers(side).y+grone(1+ismember('shift',command.Modifier))*trajvector;
%
%         end
%
%     case 'rightarrow'
%         % side=2
%         try
%             trajvector=gettraj(trajectory,2);
%             side=2;
%             newmarkers(side).head=oldmarkers(side).head-grone(1+ismember('shift',command.Modifier))*trajvector;
%             newmarkers(side).tail=oldmarkers(side).tail-grone(1+ismember('shift',command.Modifier))*trajvector;
%             newmarkers(side).x=oldmarkers(side).x-grone(1+ismember('shift',command.Modifier))*trajvector;
%             newmarkers(side).y=oldmarkers(side).y-grone(1+ismember('shift',command.Modifier))*trajvector;
%
%         end

end

switch command.Character
    case {'+','*'}
        try
            trajvector=gettraj(trajectory,options.elside);
            newmarkers(options.elside).tail=oldmarkers(options.elside).tail+trajvector.*(0.3*(ismember('shift',command.Modifier)+1));
        end

    case {'-','_'}
        try
            trajvector=gettraj(trajectory,options.elside);
            newmarkers(options.elside).tail=oldmarkers(options.elside).tail-trajvector.*(0.3*(ismember('shift',command.Modifier)+1));
        end

    case 's'
        %newmarkers(side).x=rotate(

end
%disp(command);


function traj=gettraj(trajectory,side)

traj=mean(diff(trajectory{side}));
