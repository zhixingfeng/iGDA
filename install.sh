# compile iGDA
make

# compile stxxl (require cmake)
cd tools/stxxl
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./local/stxxl ..
make install

