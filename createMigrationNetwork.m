function net = createMigrationNetwork()
% Create a network to represent the migration sample data.
% ((SPRETUS2,(SPRETUS1)I2#H1::0.6278519129107153)I6,((BACKGROUND,VARIABLE)I1,I2#H1::0.3721480870892847)I3)I0;
% ((SPRETUS2:0.9104578515297839,(SPRETUS1:2.0195704679518434)I2#H1:0.16283989808281535::0.6278519129107153)I6:2.0779118863817883,((BACKGROUND:1.146879871230814,VARIABLE:2.0600082939688695)I1:2.707813370597956,I2#H1:1.2805138209067635::0.3721480870892847)I3:0.509260077469016)I0;

% 2.0195704679518434
% 0.9104578515297839
% 0.16283989808281535
% 1.146879871230814
% 2.0600082939688695
% 2.707813370597956
% 1.2805138209067635
% 2.0779118863817883
% 0.509260077469016

buffer = calllib('libnetworkprob', 'allocNetworkBuffer', 9);
BACKGROUND = calllib('libnetworkprob', 'createLeafNetNode', buffer, 'BACKGROUND');
VARIABLE = calllib('libnetworkprob', 'createLeafNetNode', buffer, 'VARIABLE');
SPRETUS1 = calllib('libnetworkprob', 'createLeafNetNode', buffer, 'SPRETUS1');
SPRETUS2 = calllib('libnetworkprob', 'createLeafNetNode', buffer, 'SPRETUS2');

I2 = calllib('libnetworkprob', 'createNetworkNetNode', buffer, 'I2',  createEdge(0, SPRETUS1), 0.6278519129107153, 9);
I6 = calllib('libnetworkprob', 'createTreeNetNode', buffer, 'I6',  createEdge(1, SPRETUS2), createEdge(2, I2, 'LEFT'));
I1 = calllib('libnetworkprob', 'createTreeNetNode', buffer, 'I1',  createEdge(3, BACKGROUND), createEdge(4, VARIABLE));
I3 = calllib('libnetworkprob', 'createTreeNetNode', buffer, 'I3',  createEdge(5, I1), createEdge(6, I2, 'RIGHT'));
I0 = calllib('libnetworkprob', 'createTreeNetNode', buffer, 'I0',  createEdge(7, I6), createEdge(8, I3));

net.rootNode = I0;
net.buffer = buffer;

net = libstruct('Network', net);

end