# Load-balance Level Coarsening (LBC)
Load-balance Level Coarsening is a DAG partitionig/scheduling 
algorithm used for making sparse matrix loops parallel. 
It can be used within code generators or libraries. For more information see 
[Sympiler documents](https://www.sympiler.com/docs/lbc/).

## Install

### Prerequisites 
* CMake
* C++ compiler (GCC, ICC, or CLang)
* METIS (optional) dependency for running the demo efficiently.

### Linux
Then install LBC by following commands:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If METIS is installed in the system path,
CMAKE will resolve the dpendency otherwise you need to set 
`CMAKE_PREFIX_PATH` to the root directory of metis, i.e., 
where the cmakelists file exists. 
For installing METIS in Ubuntu you can also use
```
sudo apt install metis
```


### Mac
Setting the C and CXX compilers to GCC and then follow the Linux 
instructions. For example:
`-DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc\@9/9.3.0_2/bin/g++-9 -DCMAKE_C_COMPILER=/usr/local/Cellar/gcc\@9/9.3.0_2/bin/gcc-9`

Alternatively, you can use CLang using `brew install llvm`. 
The default clang on Mac might not work so make sure to set it llvm clang:
`-DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++`

## Example
The example directory shows how to call LBC API and iterate over 
the created partitioning. For more examples on how LBC is used for
making loops with sparse dependencies parallel, please check 
[Sympiler Git repo](https://github.com/sympiler/sympiler).

