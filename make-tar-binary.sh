#!/bin/sh

cd ..

mv synthlisa/examples/dev ./synthlisa-dev
mkdir ./synthlisa-data
mv synthlisa/examples/data/* ./synthlisa-data/.

tar zcvpf synthlisa/synthlisa-binary-$1-1.0.tar.gz \
    synthlisa/README.txt \
    synthlisa/Makefile \
    synthlisa/examples \
    synthlisa/setdir.sh \
    synthlisa/setdir.csh \
    synthlisa/$1

mv synthlisa-data/* synthlisa/examples/data/.
rmdir synthlisa-data
mv synthlisa-dev synthlisa/examples/dev

