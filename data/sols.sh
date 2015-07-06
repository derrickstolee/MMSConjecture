#!/bin/bash
tail -n 100000 out* | grep -f solaccept.txt 
