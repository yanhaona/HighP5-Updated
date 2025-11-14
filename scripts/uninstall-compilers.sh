#!/bin/bash

# keep track of the current directory
installer_dir=`pwd`

# list of compiler directories
multicore_compiler_dir=compilers/multicore-backend  
segmented_memory_compiler_dir=compilers/segmented-backend
frontend_compiler_dir=compilers/frontend

# common library directory
common_libs_dir=compilers/common-libs


# enter the multicore compiler directory and clean it
echo "cleaning up the IT multicore backend compiler"
cd $multicore_compiler_dir
make -f MakeFile-Compiler clean
rm -f micc
rm -f config/deployment.properties

# come back to the installer directory and delete the compiler script
cd $installer_dir
rm -f micc

# enter the segmented-memory compiler directory and clean it
echo "cleaning up the IT segmented memory backend compiler"
cd $segmented_memory_compiler_dir
make -f MakeFile-Compiler clean
rm -f sicc
rm -f config/deployment.properties

# enter the frontend compiler directory and clean it
echo "cleaning up the IT frontend compiler"
cd $installer_dir
cd $frontend_compiler_dir
make clean

# clean all common library objects
echo "cleaning up the common libraries"
cd $installer_dir
cd $common_libs_dir
find . -name "*.o" -type f -print -delete

# come back to the installer directory and delete the compiler script
cd $installer_dir
rm -f smicc
