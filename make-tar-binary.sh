#!/bin/sh

cd ..

mv synthlisa/examples/dev ./synthlisa-dev
mkdir ./synthlisa-data
mv synthlisa/examples/data/* ./synthlisa-data/.
mkdir ./synthlisa-eps
mv synthlisa/examples/eps/* ./synthlisa-eps/.
mv synthlisa/examples/CVS ./synthlisa-ecvs
mv synthlisa/examples/tdibinary-X-montana.bin .

tar zcvpf synthlisa/synthlisa-binary-$1-1.0.tar.gz \
    synthlisa/README.txt \
    synthlisa/Makefile \
    synthlisa/examples \
    synthlisa/setdir.sh \
    synthlisa/setdir.csh \
    synthlisa/$1

mv tdibinary-X-montana.bin synthlisa/examples/.
mv synthlisa-ecvs synthlisa/examples/CVS
mv synthlisa-eps/* synthlisa/examples/eps/.
rmdir ./synthlisa-eps
mv synthlisa-data/* synthlisa/examples/data/.
rmdir synthlisa-data
mv synthlisa-dev synthlisa/examples/dev
