#!/bin/sh

function usage {
    echo "Usage: $0 <destination base directory>"
    echo "The destination base directory must be writable"
    exit 0
}

if [ "x$1" == "x" ] || ! [ -d $1 ]; then
    usage
fi

CHROMOSOMES="01 02 03 04 05 06 07 08 09 10 11 12"
FINDIRS="finished unfinished"
SUBDIRS="annotation/";
for chrnum in $CHROMOSOMES; do
    for fin in $FINDIRS; do
	for subdir in $SUBDIRS; do
	    mkdir -p $1/chr$chrnum/$fin/$subdir
	done
    done
done
