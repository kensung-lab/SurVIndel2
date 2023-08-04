tar -xjf htslib-1.18.tar.bz2
cd htslib-1.18
autoheader
autoconf
./configure --prefix=`pwd`
make
make install
