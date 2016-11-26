function edge = createEdge(index, sourceNode, type)
% Create an edge with the corresponding edge index, source node and type.
% The type is optional and defaults to a normal edge.
% See the corresponding C++ code for better documentation of what the edge types mean.

    if nargin == 2
        type = 'NORMAL';
    end

    result.index = index;
    result.sourceNode = sourceNode;
    result.distance = 0;
    result.type = type;

    edge = libstruct('NetworkEdge', result);
end