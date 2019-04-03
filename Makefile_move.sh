#!/bin/sh

cuda_capability=$(deviceQuery | grep CUDA\ Capability | cut -d " " -f 11)
if [ $cuda_capability = 2.0 ]; then
	echo 2.0
	#cp Makefile.cc20 Makefile
elif [ $cuda_capability = 6.1 ]; then
	echo 6.1
	#cp Makefile.cc61 Makefile
elif [ $cuda_capability = 7.5 ]; then
	echo 7.5
	#cp Makefile.cc75 Makefile
else
	exit 1
fi


