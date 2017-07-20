c -c prepare.cu --compiler-options -fPIC
gcc -o debug.o -c debug.cpp
gcc -o debug debug.o prepare.o -lcudart -L"/usr/local/cuda-8.0/lib64" -lm
