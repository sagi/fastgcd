#!/bin/bash

__exists() {
    which $1 1>/dev/null 2>&1
}

get="fetch";
! __exists fetch && get="curl -OL";

if [ ! -d gmp-5.0.5 ];
then

    if [ ! -f gmp-5.0.5.tar.bz2 ];
    then
        $get ftp://ftp.gmplib.org/pub/gmp-5.0.5/gmp-5.0.5.tar.bz2
    fi

    sum=`openssl sha1 gmp-5.0.5.tar.bz2 | awk -F' ' '{print $2}'`

    if [[ $sum != "12a662456033e21aed3e318aef4177f4000afe3b" ]];
    then
        echo ''
        echo '=========================================='
        echo 'ERROR: could not verify gmp-5.0.5.tar.bz2;'
        echo 'Downloaded over untrusted channel (non-https)'
        echo '=========================================='
        exit;
    fi

    tar xf gmp-5.0.5.tar.bz2
fi

cd gmp-5.0.5
patch -p 1 < ../gmp-5.0.5.patch
mkdir ../gmp-patched
./configure --prefix=$PWD/../gmp-patched/
make
make install
cd ..
make

