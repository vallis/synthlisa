#!/bin/sh

cd ..

mv synthlisa/examples/dev ./synthlisa-dev
mkdir ./synthlisa-data
mv synthlisa/examples/data/* ./synthlisa-data/.
mkdir ./synthlisa-eps
mv synthlisa/examples/eps/* ./synthlisa-eps/.
mv synthlisa/examples/CVS ./synthlisa-ecvs
#mv synthlisa/examples/tdibinary-X-montana.bin .
mv synthlisa/lisasim/build ./synthlisa-build
mv synthlisa/lisasim/dev ./synthlisa-ldev
mv synthlisa/lisasim/CVS ./synthlisa-CVS
mv synthlisa/contrib-source/GSL-1.4/build ./synthlisa-gbuild
mv synthlisa/contrib-source/GSL-1.4/CVS ./synthlisa-gcvs

rm synthlisa/examples/.DS_Store
rm synthlisa/lisasim/.DS_Store

tar zcvpf synthlisa/synthlisa-source-1.0.tar.gz \
    synthlisa/README.txt \
    synthlisa/Makefile \
    synthlisa/examples \
    synthlisa/setdir.sh \
    synthlisa/setdir.csh \
    synthlisa/make-tar-binary.sh \
    synthlisa/lisasim \
    synthlisa/contrib-source/Makefile \
    synthlisa/contrib-source/GSL-1.4 \
    synthlisa/contrib-source/tars \
    synthlisa/license.pdf \
    synthlisa/license.doc

mv synthlisa-gcvs synthlisa/contrib-source/GSL-1.4/CVS
mv synthlisa-gbuild synthlisa/contrib-source/GSL-1.4/build
mv synthlisa-CVS synthlisa/lisasim/CVS
mv synthlisa-ldev synthlisa/lisasim/dev
mv synthlisa-build synthlisa/lisasim/build
#mv tdibinary-X-montana.bin synthlisa/examples/.
mv synthlisa-ecvs synthlisa/examples/CVS
mv synthlisa-eps/* synthlisa/examples/eps/.
rmdir ./synthlisa-eps
mv synthlisa-data/* synthlisa/examples/data/.
rmdir synthlisa-data
mv synthlisa-dev synthlisa/examples/dev
