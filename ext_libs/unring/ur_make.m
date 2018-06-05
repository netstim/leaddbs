if isunix
    compflags = ' COMPFLAGS=''$COMPFLAGS -static-libstdc++'' CXXFLAGS=''$CXXFLAGS -msse2''';
    fftwlib = 'libfftw3.a';
elseif ispc
    compflags = ' COMPFLAGS="$COMPFLAGS /MT" CXXFLAGS="$CXXFLAGS /arch:AVX"';
    fftwlib = 'fftw3.lib';
end

eval(['mex -compatibleArrayDims' compflags, ' ringRm.cpp ', fftwlib]);
