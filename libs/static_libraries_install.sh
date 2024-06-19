#!/bin/bash
set -e

# apt install libboost-all-dev

wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
tar -xzf boost_1_81_0.tar.gz
cd boost_1_81_0
./bootstrap.sh
./b2 --build-dir=build_temp link=static --with-random
./b2 install --prefix=../Boost --with-random link=static
cd ..
rm -R boost_1_81_0

echo ""

sudo apt install cmake

echo ""

git clone https://github.com/BioPP/bpp-core

cd bpp-core
mkdir build
cd build
cmake -B. -S.. -DCMAKE_INSTALL_PREFIX=../../Bpp -DBUILD_STATIC=ON  # prepare compilation
make # compile
make install # move files to the installation directory

cd ../..

rm -R bpp-core

echo ""

git clone https://github.com/BioPP/bpp-seq

cd bpp-seq
mkdir build
cd build
cmake -B. -S.. -DCMAKE_INSTALL_PREFIX=../../Bpp -DBUILD_STATIC=ON # prepare compilation
make # compile
make install # move files to the installation directory

cd ../..
rm -R bpp-seq

echo ""

apt install libeigen3-dev

echo ""

git clone https://github.com/BioPP/bpp-phyl

cd bpp-phyl
mkdir build
cd build
cmake -B. -S.. -DCMAKE_INSTALL_PREFIX=../../Bpp  -DBUILD_STATIC=ON # prepare compilation
make # compile
make install # move files to the installation directory

cd ../..
rm -R bpp-phyl
