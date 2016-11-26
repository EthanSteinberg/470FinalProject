function trees = createAllTrees(leafNodes)
% Create all possible binary trees from a given list of leaf nodes.


children = createParentMapForTrees(length(leafNodes));
totalNodes = length(children) * (length(leafNodes) * 2 - 1);

buffer = calllib('libnetworkprob', 'allocTreeBuffer', totalNodes);

for i = 1:size(children, 1)
    root = children(i, 1, 1);
    trees(i).rootNode = createTree(buffer, root, squeeze(children(i, :, :)), leafNodes);
    trees(i).buffer = buffer;
end

end

function tree = createTree(buffer, node, children, leafNodes)
    if node > 0
        tree = calllib('libnetworkprob', 'createLeafTreeNode', buffer, leafNodes(node));
    elseif node < 0
        leftChild = createTree(buffer, children(-node, 1), children, leafNodes);
        rightChild = createTree(buffer, children(-node, 2), children, leafNodes);
        tree = calllib('libnetworkprob', 'createTreeTreeNode', buffer, num2str(-node), leftChild, rightChild);
    end
end


function children = createParentMapForTrees(numberOfLeaves)
edges = {{{1, -1}}};
for i = 2:numberOfLeaves
    nextEdges = {};

    for treeIndex = 1:length(edges)
       tree = edges{treeIndex};

       for edge = 1:length(tree)
           newEdges = tree;
           newEdges{edge} = {tree{edge}{1}, -i};
           newEdges{length(newEdges) + 1} = {-i, tree{edge}{2}};
           newEdges{length(newEdges) + 1} = {i, -i};

           nextEdges{length(nextEdges) + 1} = newEdges;
       end
    end

    edges = nextEdges;
end

for treeIndex = 1:length(edges)

    for i = 1:numberOfLeaves
        children(treeIndex, i, :) = [0, 0];
    end

    tree = edges{treeIndex};
    for edgeIndex = 1:length(tree)
        source = tree{edgeIndex}{1};
        destination = tree{edgeIndex}{2};

        if children(treeIndex, -destination, 1) == 0
            children(treeIndex, -destination, 1) = source;
        else
            children(treeIndex, -destination, 2) = source;
        end
    end
end

end




