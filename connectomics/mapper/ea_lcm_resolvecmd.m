function cmds=ea_lcm_resolvecmd(cmd)
switch cmd
    case 1
        cmds='seed';
    case 2
        cmds='pseed';
    case 3
        cmds='matrix';
    case 4
        cmds='pmatrix';
    case 5
        cmds='pmap';
end
