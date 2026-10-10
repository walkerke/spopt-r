# Contributing to spopt

Thank you for your interest in contributing to spopt! This document
provides instructions for setting up a development environment and
building the package from source.

## Installation for users

Most users should install pre-built binaries from r-universe:

\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"spopt"``, repos ``=`` ``"https://walkerke.r-universe.dev"``)`

The instructions below are for **contributors and developers** who need
to build from source.

## Building from source

spopt uses a Rust backend via [extendr](https://extendr.github.io/) and
the [HiGHS](https://highs.dev/) solver for mixed-integer programming.
Building from source requires several system dependencies.

### macOS

1.  **Install Rust** (if not already installed):

    ``` bash
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    source $HOME/.cargo/env
    ```

2.  **Install system dependencies** via Homebrew:

    ``` bash
    brew install cmake llvm
    ```

3.  **Install the package**:

    \
    `# install.packages("pak")`\
    `pak``::`[`pak`](https://pak.r-lib.org/reference/pak.html)`(``"walkerke/spopt-r"``)`

### Linux (Ubuntu/Debian)

1.  **Install Rust**:

    ``` bash
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    source $HOME/.cargo/env
    ```

2.  **Install system dependencies**:

    ``` bash
    sudo apt-get update
    sudo apt-get install -y cmake libclang-dev
    ```

3.  **Install the package**:

    \
    `pak``::`[`pak`](https://pak.r-lib.org/reference/pak.html)`(``"walkerke/spopt-r"``)`

### Windows

Windows builds are more complex due to toolchain requirements. We
strongly recommend using r-universe binaries on Windows.

If you must build from source:

1.  **Install Rtools44** (or appropriate version for your R):

    - Download from <https://cran.r-project.org/bin/windows/Rtools/>

2.  **Install Rust** using rustup:

    - Download from <https://rustup.rs/>
    - During installation, select the `x86_64-pc-windows-gnu` toolchain

3.  **Configure Rust for GNU toolchain**:

    ``` bash
    rustup default stable-x86_64-pc-windows-gnu
    ```

4.  **Install additional MSYS2 packages** (from Rtools terminal):

    ``` bash
    pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-clang
    ```

5.  **Install the package**:

    \
    `pak``::`[`pak`](https://pak.r-lib.org/reference/pak.html)`(``"walkerke/spopt-r"``)`

## Development workflow

### Making changes to R code

After modifying R code in `R/`, regenerate documentation:

\
`devtools``::`[`document`](https://devtools.r-lib.org/reference/document.html)`(``)`\
`devtools``::`[`check`](https://devtools.r-lib.org/reference/check.html)`(``)`

### Making changes to Rust code

The Rust source is in `src/rust/`. After making changes:

\
`# Recompile Rust code and reload`\
`rextendr``::`[`document`](https://extendr.github.io/rextendr/reference/document.html)`(``)`\
\
`# Or for a full rebuild`\
`devtools``::`[`load_all`](https://devtools.r-lib.org/reference/load_all.html)`(``)`

### Running tests

\
`devtools``::`[`test`](https://devtools.r-lib.org/reference/test.html)`(``)`

### Building the pkgdown site

The vignettes use Quarto. To build the documentation site:

\
`pkgdown``::`[`build_site`](https://pkgdown.r-lib.org/reference/build_site.html)`(``)`

## Reporting issues

Please report bugs and feature requests at
<https://github.com/walkerke/spopt-r/issues>.

When reporting bugs, please include: - Your operating system and
version - R version
([`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)) - A
minimal reproducible example
