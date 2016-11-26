function net = createSimpleNetwork()
% Create a simple network just for testing purposes.
% It consists of three leaves and one introgression point in the middle.

buffer = calllib('libnetworkprob', 'allocNetworkBuffer', 7);
A = calllib('libnetworkprob', 'createLeafNetNode', buffer, 'A');
B = calllib('libnetworkprob', 'createLeafNetNode', buffer, 'B');
C = calllib('libnetworkprob', 'createLeafNetNode', buffer, 'C');

BIntrogressed =  calllib('libnetworkprob', 'createNetworkNetNode', buffer, 'C',  createEdge(0, B), 0.5, 7);

one = calllib('libnetworkprob', 'createTreeNetNode', buffer, 'one',  createEdge(1, A), createEdge(2, BIntrogressed, 'LEFT'));
two = calllib('libnetworkprob', 'createTreeNetNode', buffer, 'two',  createEdge(3, BIntrogressed, 'RIGHT'), createEdge(4, C));
three = calllib('libnetworkprob', 'createTreeNetNode', buffer, 'three',  createEdge(5, one), createEdge(6, two));

net.rootNode = three;
net.buffer = buffer;

net = libstruct('Network', net);

end