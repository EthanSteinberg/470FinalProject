function tree = generateTreeFromNewick(newick)
% Create a tree from a given newick string.

numLeafs = sum(newick == ',') + 1;

totalNodes = numLeafs * 2 - 1;

tree.buffer = calllib('libnetworkprob', 'allocTreeBuffer', totalNodes);

tree.rootNode = getNode(tree.buffer, newick, 1);

end

function [node, finalEnd] = getNode(buffer, newick, startIndex)

if newick(startIndex) == '('
    % This is a compound node.
    [leftChild, middle] = getNode(buffer, newick, startIndex+1);
    [rightChild, finalEnd] = getNode(buffer, newick, middle+2);
    finalEnd = finalEnd + 1;
    node = calllib('libnetworkprob', 'createTreeTreeNode', buffer, num2str(startIndex), leftChild, rightChild);
else
    finalEnd = find(newick(startIndex:end) == ',' | newick(startIndex:end) == ')',1,'first') + startIndex - 2;
    node = calllib('libnetworkprob', 'createLeafTreeNode', buffer, newick(startIndex:finalEnd));
end
end


