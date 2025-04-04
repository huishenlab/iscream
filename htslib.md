## Getting htslib headers

You or your system administrator can install htslib with the system package
manager which usually sets `PKG_CONFIG_PATH` automatically. On MacOS you can get
htslib with the Homebrew package manager. On HPC systems, htslib may be provided
as a module. Make sure these methods also set the `PKG_CONFIG_PATH`.

If you aren't able to install htslib development libraries system-wide for lack
of admin permissions, you can install them from other channels. We recommend
[pixi](#pixi) or [conda](#conda) since their htslib is compiled with
[*libdeflate*](https://github.com/ebiggers/libdeflate) support and is faster
than htslib without *libdeflate*. If you're compiling your own htslib, compile
*libdeflate* first and then htslib. You can also use Rhtslib, but we've found
this has poorer performance than a standard htslib installation, both with and
without *libdeflate*.

### Conda/miniconda/mamba/micromamba

The *iscream* repo provides [environment.yaml](environment.yaml) to install
htslib 1.21 - add this file to your project directory and run

```bash
conda env create -f environment.yaml
conda activate iscream
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib/
```

Confirm that the headers are available for compilation
```bash
pkg-config --cflags --libs htslib
```

You should get something like

```
-I/home/user/miniconda3/envs/iscream/include -L/home/user/miniconda3/envs/iscream/lib -lhts
```

### [Pixi](https://pixi.sh/latest/)

Pixi uses the conda repositories to install packages. The *iscream* repo
provides [pixi.toml](pixi.toml) to install htslib 1.21 - add this file to your
project directory and run

```bash
pixi shell
```

Confirm that the headers are available for compilation
```bash
pkg-config --cflags --libs htslib
```

You should get something like

```
-I/home/user/iscream/.pixi/envs/default/include -L/home/user/iscream/.pixi/envs/default/lib -lhts
```

### Rhtslib

If the brew/conda/pixi method doesn't work, you can use the Rhtslib branch to
get the headers. As Rats lib's headers are a few years older than the latest
htslib version, this branch may not be consistently supported as iscream updates
to follow the latest htslib API.

```
devtools::install_github("huishenlab/iscream@rhtslib")
pak::pkg_install("huishenlab/iscream@rhtslib")
```
