tar -xjf htslib-1.21.tar.bz2
cd htslib-1.21
autoheader
autoconf
./configure --prefix=`pwd`
make
make install
