#!/bin/sh

cd ..

mv synthlisa/examples/dev ./synthlisa-dev
mkdir ./synthlisa-data
mv synthlisa/examples/data/* ./synthlisa-data/.
mv synthlisa/lisasim/build ./synthlisa-build
mv synthlisa/lisasim/dev ./synthlisa-ldev
mv synthlisa/lisasim/CVS ./synthlisa-CVS

tar zcvpf synthlisa/synthlisa-source-1.0.tar.gz \
    synthlisa/README.txt \
    synthlisa/Makefile \
    synthlisa/examples \
    synthlisa/setdir.sh \
    synthlisa/setdir.csh \
    synthlisa/lisasim \
    synthlisa/contrib-source/Makefile \
    synthlisa/contrib-source/GSL-1.4 \
    synthlisa/contrib-source/tars

mv synthlisa-CVS synthlisa/lisasim/CVS
mv synthlisa-ldev synthlisa/lisasim/dev
mv synthlisa-build synthlisa/lisasim/build
mv synthlisa-data/* synthlisa/examples/data/.
rmdir synthlisa-data
mv synthlisa-dev synthlisa/examples/dev
