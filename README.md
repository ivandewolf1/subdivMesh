# graded subdivision meshes

this document will soon be written in full.
If nothing else, jump to the bottom and click on the youtube link.

## compiling
***These instructions are for linux***

### step 1
this uses the lightwight templated vector library, you need the 3 .h files

https://github.com/ivandewolf1/lwtv

I made a sympolic link to the .h files:
```
ln -s /job/git/lwtv/src/tvec3.h  src/lwtv/tvec3.h 
ln -s /job/git/lwtv/src/tmat3.h  src/lwtv/tmat3.h 
ln -s /job/git/lwtv/src/tmat4.h  src/lwtv/tmat4.h
``` 
but, it may be easist to simply copy those 3 files in.

### step 2
#### setup Houdini
this has a sample SOP for Houdini as it's testbed. 
A free version of Houdini is available here:
https://www.sidefx.com/

this is currently being tested in Houdini 17.5.229
```
pushd /opt/hfs17.5
source houdini_setup
popd
```

### step 3
#### cmake
create a build directory, then create a makefile with cmake
```
mkdir build
cd build
cmake ../src
```

### step 4
#### make
assuming that you are on linux, and cmake ran correctly, you will probably be set up now to run make.
```
make
```
### testing
if it compiled correctly, it will have made a .so file in your home directory.
```
prompt$ ls ~/houdini17.5/dso
SOP_SubtriBasic.so
```
this youtube video will show you how to see if the node loads into Houdini:

https://youtu.be/Sw7N1YM6Dnc
