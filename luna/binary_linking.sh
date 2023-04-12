install_name_tool -add_rpath @executable_path/. luna
install_name_tool -change /opt/homebrew/opt/fftw/lib/libfftw3.3.dylib @rpath/libfftw3.3.dylib luna
install_name_tool -change /opt/homebrew/opt/libomp/lib/libomp.dylib @rpath/libomp.dylib luna
install_name_tool -change /opt/homebrew/opt/lightgbm/lib/lib_lightgbm.dylib @rpath/lib_lightgbm.dylib luna
