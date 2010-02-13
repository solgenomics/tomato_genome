#!/bin/sh
ftpdir='/data/shared/ftpsite/tomato_genome/bacs'
uploaddir='/data/shared/tomato_genome/country_uploads/'

#pattern that excludes any upload directories
#that we don't process for submissions
exclude_pattern='/fishstack/|/eileen/|/fishchina/|/SKEL/'

processed_dirs=`find $uploaddir -type d -and -name upload-processed | egrep -v $exclude_pattern`;

for dir in $processed_dirs; do
    cd $dir;
    sudo mv -v *.tar.gz ../upload/;
    cd -;
done
