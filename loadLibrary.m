function loadLibrary()
% A simple utility for loading my network probability code

% Point to the directory where the code is located
directory = 'networkprob';

if libisloaded('libnetworkprob')
    unloadlibrary('libnetworkprob');
end

loadlibrary([directory '/build/libnetworkprob'], [directory '/src/matlabffi.h']);
