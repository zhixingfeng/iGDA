# compile stxxl (require cmake)
echo "================ install iGDA dependences ==============="
echo "install stxxl"
cd tools/stxxl
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./local/stxxl .. > build.log
make install >> build.log
cd ../../..

# configure stxxl
echo "disk=./tmp_stxxl,32G,syscall unlink" > $HOME/.stxxl

# compile graphviz
echo "install graphviz"
cd tools/graphviz/graphviz
mkdir -p build
./configure --prefix=$(pwd)/build > build.log
make >> build.log
make install >> build.log

cd ../../..

# compile iGDA
make

# copy tred into bin folder
cp ./tools/graphviz/graphviz/build/bin/tred bin/


