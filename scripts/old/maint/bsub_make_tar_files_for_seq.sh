#!/bin/sh

for seqfile in $@; do
    basename=`basename $seqfile .seq`
    mkdir /tmp/$basename;
    cp $seqfile /tmp/$basename/;
    tar -C /tmp -czf $basename.tar.gz $basename/
    rm -rf /tmp/$basename/;
done