#!/usr/bin/csh

setenv fileName $argv[1]


setenv nline `wc -l $fileName | awk '{print($1)}'`

awk -F: -v nl="$nline" '(NR<nl) {print($0)}' $fileName > tmp
mv tmp $fileName







