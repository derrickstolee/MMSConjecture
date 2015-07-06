#!/bin/bash

../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N  5 -K 4 < in.0 > out-04-05-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N  6 -K 4 < in.0 > out-04-06-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N  7 -K 4 < in.0 > out-04-07-$1-$2.txt
#../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N  8 -K 4 < in.0 > out-04-08-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N  9 -K 4 < in.0 > out-04-09-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 10 -K 4 < in.0 > out-04-10-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 11 -K 4 < in.0 > out-04-11-$1-$2.txt
#../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 12 -K 4 < in.0 > out-04-12-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 13 -K 4 < in.0 > out-04-13-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 14 -K 4 < in.0 > out-04-14-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 15 -K 4 < in.0 > out-04-15-$1-$2.txt
#../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 16 -K 4 < in.0 > out-04-16-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 17 -K 4 < in.0 > out-04-17-$1-$2.txt
