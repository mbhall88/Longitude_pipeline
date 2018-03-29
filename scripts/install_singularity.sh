#!/usr/bin/env sh
VERSION=2.4.5
wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-${VERSION}.tar.gz
tar xvf singularity-${VERSION}.tar.gz
cd singularity-${VERSION}
./configure --prefix=/usr/local
make
make install
rm singularity-${VERSION}.tar.gz
