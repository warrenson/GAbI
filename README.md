# GAbI
Genome Annotation by Information

## Environment

### Build
Requires [pixi](https://pixi.prefix.dev/latest/installation/)

```
pixi install
```

### Fast BiEntropy
Need to build from source to get the C extension impelmentation (faster).

Here's how I did it, YMMV.

1. Requires `libgmp` provided by `gmp` in homebrew
2. Export the `gmp` lib path for clang on MacOS (M1) 13.5.2

```
export LIBRARY_PATH=/opt/homebrew/Cellar/gmp/6.3.0/lib
export CPATH=/opt/homebrew/Cellar/gmp/6.3.0/include
```

3. Get the source, and install it via the pixi env

```
curl -sSL https://github.com/sandialabs/bientropy/archive/refs/tags/v1.1.4.tar.gz | tar xz
cd bientropy-1.1.4
pixi run --manifest-path $(readlink -f ..) python setup.py install
```
