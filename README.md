## Installation

- Compile the code using cmake>3.20.0

```bash
cd wildmeshing-toolkit
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

You may need to install `gmp` before compiling the code. You can install `gmp` via [homebrew](https://brew.sh/).

## Usage

To reproduce figures from the paper, please use the commands from [reproduce_scripts](reproduce_scripts.sh).

If the script does not execute properly, it might be because of CRLF line endings. You can fix that by running `dos2unix reproduce_scripts.sh`.