set -e
sudo -E apt update 
sudo -E apt-get install -yq --no-install-recommends build-essential \
 wget cmake  libboost-filesystem-dev libboost-program-options-dev libboost-iostreams-dev libboost-date-time-dev \
 libprotoc-dev libprotoc-dev protobuf-compiler \
 mafft rsync libtbb-dev openmpi-bin libopenmpi-dev automake libtool autoconf make nasm ninja-build

# create build directory
startDir=$pwd
cd $(dirname "$0")
mkdir -p ../build
cd ../build
wget https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz
tar -xvf v2.30.0.tar.gz
cd isa-l-2.30.0
./autogen.sh
./configure
make -j$(nproc)
sudo -E make install
cd ..

#download and install TBB
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz 
tar -xvzf 2019_U9.tar.gz

#download and install fmt library
FMT_VERSION=8.1.1
curl -OL https://github.com/fmtlib/fmt/archive/$FMT_VERSION.tar.gz
tar -xvf $FMT_VERSION.tar.gz

cmake -GNinja -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  -DFMT_DIR=${PWD}/fmt-$FMT_VERSION ..
#make -j$(nproc) VERBOSE=1
ninja

# install faToVcf
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
chmod +x faToVcf

cd $startDir