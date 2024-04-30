function tree = ea_struct2tree(S, tree, options)
% Construct uitree based on input struct

arguments
    S struct
    tree matlab.ui.container.Tree = uitree
    options.expandAll logical {mustBeNumericOrLogical} = false
end

cellfun(@(f) addNode(tree, S, f), fieldnames(S));

if options.expandAll
    expand(tree, 'all');
end


function addNode(parentNode, S, field)
childNode = uitreenode(parentNode, 'Text', field);
if isstruct(S.(field))
    cellfun(@(f) addNode(childNode, S.(field), f), fieldnames(S.(field)));
else
    childNode.NodeData = S.(field);
end
