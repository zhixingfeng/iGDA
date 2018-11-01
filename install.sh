# compile stxxl (require cmake)
echo "================ install iGDA dependences ==============="
echo "install stxxl"
cd tools/stxxl
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./local/stxxl ..
make install
cd ../../..

# configure stxxl
echo "disk=./tmp_stxxl,32G,syscall unlink" > $HOME/.stxxl

# compile graphviz
echo "install graphviz"
cd tools/graphviz/graphviz-2.40.1
./autogen.sh
mkdir -p build
cd build
../configure --prefix=$(pwd)
make
make install

cd ../../../..

# compile iGDA
echo "build iGDA"
make

# copy tred into bin folder
cp ./tools/graphviz/graphviz-2.40.1/build/bin/tred bin/igda_tred


