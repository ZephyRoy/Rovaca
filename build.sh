mkdir build
mkdir release
cd build
cmake -DCMAKE_PREFIX_PATH=/usr/include/boost169/ -DBOOST_ROOT=/usr/include/boost169/ -DCMAKE_INSTALL_PREFIX=../release ..
make -j90 && make install
cd ..
strip -s release/lib/libIOStream.so  release/lib/libassemble.so  release/lib/libbase.so  release/lib/libgenotype.so  release/lib/libhaplotypecaller.so  release/lib/liblogger.so  release/lib/libpairhmm.so  release/lib/libregion.so  release/lib/libutils.so  release/lib/libwriter.so release/lib/libbqsr.so release/lib/libmain.so release/lib/rovaca