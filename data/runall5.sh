#!/bin/bash

../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N  6 -K 5 < in.0 > out-05-06-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N  7 -K 5 < in.0 > out-05-07-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N  8 -K 5 < in.0 > out-05-08-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N  9 -K 5 < in.0 > out-05-09-$1-$2.txt
#../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 10 -K 5 < in.0 > out-05-10-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 11 -K 5 < in.0 > out-05-11-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 12 -K 5 < in.0 > out-05-12-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 13 -K 5 < in.0 > out-05-13-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 14 -K 5 < in.0 > out-05-14-$1-$2.txt
#../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 15 -K 5 < in.0 > out-05-15-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 16 -K 5 < in.0 > out-05-16-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 17 -K 5 < in.0 > out-05-17-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 18 -K 5 < in.0 > out-05-18-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 19 -K 5 < in.0 > out-05-19-$1-$2.txt
#../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 20 -K 5 < in.0 > out-05-20-$1-$2.txt
../src/mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --rule $2 -N 21 -K 5 < in.0 > out-05-21-$1-$2.txt 

