Complie STRUMPACK on Linux
I use Centos8 


- Install metis

wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
tar -xvzf metis-5.1.0.tar.gz
cd metis-5.1.0
mkdir install
make config shared=1 cc=gcc prefix=`pwd`/install       # specify your C compiler as cc=..
make
make install
export METIS_DIR=`pwd`/install    # this will help STRUMPACK find metis


- Get ZFP
git clone https://github.com/LLNL/zfp.git
cd zfp
mkdir build
mkdir install
cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install
make -j8
make install
cd ../
export ZFP_DIR=`pwd`/install



git clone https://github.com/pghysels/STRUMPACK.git
cd STRUMPACK
git checkout matlab

mkdir build
mkdir install
cd build

cmake ../ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=../install \
  -DSTRUMPACK_USE_MPI=OFF \
  -DTPL_BLAS_LIBRARIES="/home/mingyu003/OpenBLAS/libopenblas.a"\
  -DTPL_LAPACK_LIBRARIES="/home/mingyu003/OpenBLAS/libopenblas.a"\
  -DMatlab_ROOT_DIR=/usr/local/MATLAB/R2019a/ \
  -DBUILD_SHARED_LIBS=ON \
  -DTPL_ENABLE_MATLAB=ON \
  -DTPL_ENABLE_ZFP=ON
make -j32
make install
