# compile stxxl (require cmake)
#echo "================ install iGDA dependences ==============="
#echo "install stxxl"
#cd tools/stxxl
#mkdir -p build
#cd build
#cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./local/stxxl ..
#make install
#cd ../../..

# configure stxxl
#echo "disk=./tmp_stxxl,8G,syscall unlink" > $HOME/.stxxl

# compile graphviz
#echo "install graphviz"
#cd tools/graphviz/
#tar -zxvf graphviz.tar.gz
#cd ./graphviz-2.40.1
#./autogen.sh
#mkdir -p build
#./configure --prefix=$(pwd)/build
#make
#make install

#cd ../../..

# compile iGDA
echo "build iGDA"
make

# copy script to bin
#cp ./script/igda_pipe_ont ./bin

#chmod u+x ./bin/igda_pipe_ont

# copy tred into bin folder
#cp ./tools/graphviz/graphviz-2.40.1/build/bin/tred bin/igda_tred


