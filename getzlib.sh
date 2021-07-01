wget --no-clobber http://sourceforge.net/projects/libpng/files/zlib/1.2.8/zlib-1.2.8.tar.gz
tar zxvf zlib-1.2.8.tar.gz
cd zlib-1.2.8
wget --no-check-certificate https://raw.githubusercontent.com/Alexpux/MSYS2-packages/master/zlib/1.2.7-minizip-cygwin.patch
wget --no-check-certificate https://raw.githubusercontent.com/Alexpux/MSYS2-packages/master/zlib/1.2.7-zlib-symbols.patch
patch -p2 -i 1.2.7-minizip-cygwin.patch
patch -p2 -i 1.2.7-zlib-symbols.patch
perl -i -wpe s/zlib1/libz/ win32/Makefile*
perl -i -wpe s/zlib1/libz/ Makefile*
perl -i -wpe s/zlib1/libz/ configure
perl -i -wpe 's/-O2/-Ofast -march=native -fPIC/' `find . -name Makefile\*`
bash configure
make
cd ..
mkdir -p z
cd z
ln -sf ../zlib-1.2.8/*.c .
ln -sf ../zlib-1.2.8/*.h .
cd ..
