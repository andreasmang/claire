# build libaries and code

cd external
./build_libs.sh --build
cd ..
make -j


# run registration

./bin/runcoldreg
