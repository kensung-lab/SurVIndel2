tar -xjf htslib-1.17.tar.bz2
cd htslib-1.17
autoheader
autoconf
./configure --prefix=`pwd`
make
make install
