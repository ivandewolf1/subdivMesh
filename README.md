# graded subdivision meshes

this document will soon be written in full.
If nothing else, jump to the bottom and click on the youtube link.

## compiling
***These instructions are for linux***

### step 1
this uses the lightwight templated vector library

https://github.com/ivandewolf1/lwtv

I cloned lwtv into the same directory that I cloned subdivMesh into, and the cmake file is expecting that.

### step 2
#### cmake
create a build directory, then create a makefile with cmake
```
mkdir build
cd build
cmake ../src
```

### step 3
#### make
assuming that you are on linux, and cmake ran correctly, you will probably be set up now to run make.
```
make
```
### testing
if it compiled correctly, it will have made an executeable in the build directory. Run this exeecuteable, and it will create an .obj geometry file
```
./tester
cat output.obj 
```

