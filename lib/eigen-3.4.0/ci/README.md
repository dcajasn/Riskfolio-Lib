## Eigen CI infrastructure

Eigen's CI infrastructure uses two stages: A `build` stage to build the unit-test
suite and a `test` stage to run the unit-tests.

### Build Stage

The build stage consists of the following jobs:

| Job Name                                 | Arch      | OS             | Compiler   | C++11   |
|------------------------------------------|-----------|----------------|------------|---------|
| `build:x86-64:linux:gcc-4.8:cxx11-off`   | `x86-64`  | `Ubuntu 18.04` | `GCC-4.8`  | `Off`   |
| `build:x86-64:linux:gcc-4.8:cxx11-on`    | `x86-64`  | `Ubuntu 18.04` | `GCC-4.8`  | `On`    |
| `build:x86-64:linux:gcc-9:cxx11-off`     | `x86-64`  | `Ubuntu 18.04` | `GCC-9`    | `Off`   |
| `build:x86-64:linux:gcc-9:cxx11-on`      | `x86-64`  | `Ubuntu 18.04` | `GCC-9`    | `On`    |
| `build:x86-64:linux:gcc-10:cxx11-off`    | `x86-64`  | `Ubuntu 18.04` | `GCC-10`   | `Off`   |
| `build:x86-64:linux:gcc-10:cxx11-on`     | `x86-64`  | `Ubuntu 18.04` | `GCC-10`   | `On`    |
| `build:x86-64:linux:clang-10:cxx11-off`  | `x86-64`  | `Ubuntu 18.04` | `Clang-10` | `Off`   |
| `build:x86-64:linux:clang-10:cxx11-on`   | `x86-64`  | `Ubuntu 18.04` | `Clang-10` | `On`    |
| `build:aarch64:linux:gcc-10:cxx11-off`   | `AArch64` | `Ubuntu 18.04` | `GCC-10`   | `Off`   |
| `build:aarch64:linux:gcc-10:cxx11-on`    | `AArch64` | `Ubuntu 18.04` | `GCC-10`   | `On`    |
| `build:aarch64:linux:clang-10:cxx11-off` | `AArch64` | `Ubuntu 18.04` | `Clang-10` | `Off`   |
| `build:aarch64:linux:clang-10:cxx11-on`  | `AArch64` | `Ubuntu 18.04` | `Clang-10` | `On`    |

### Test stage

In principle every build-job has a corresponding test-job, however testing supported and unsupported modules is divided into separate jobs. The test jobs in detail:

### Job dependecies

| Job Name                                            | Arch      | OS             | Compiler   | C++11   | Module
|-----------------------------------------------------|-----------|----------------|------------|---------|--------
| `test:x86-64:linux:gcc-4.8:cxx11-off:official`      | `x86-64`  | `Ubuntu 18.04` | `GCC-4.8`  | `Off`   | `Official`
| `test:x86-64:linux:gcc-4.8:cxx11-off:unsupported`   | `x86-64`  | `Ubuntu 18.04` | `GCC-4.8`  | `Off`   | `Unsupported`
| `test:x86-64:linux:gcc-4.8:cxx11-on:official`       | `x86-64`  | `Ubuntu 18.04` | `GCC-4.8`  | `On`    | `Official`
| `test:x86-64:linux:gcc-4.8:cxx11-on:unsupported`    | `x86-64`  | `Ubuntu 18.04` | `GCC-4.8`  | `On`    | `Unsupported`
| `test:x86-64:linux:gcc-9:cxx11-off:official`        | `x86-64`  | `Ubuntu 18.04` | `GCC-9`    | `Off`   | `Official`
| `test:x86-64:linux:gcc-9:cxx11-off:unsupported`     | `x86-64`  | `Ubuntu 18.04` | `GCC-9`    | `Off`   | `Unsupported`
| `test:x86-64:linux:gcc-9:cxx11-on:official`         | `x86-64`  | `Ubuntu 18.04` | `GCC-9`    | `On`    | `Official`
| `test:x86-64:linux:gcc-9:cxx11-on:unsupported`      | `x86-64`  | `Ubuntu 18.04` | `GCC-9`    | `On`    | `Unsupported`
| `test:x86-64:linux:gcc-10:cxx11-off:official`       | `x86-64`  | `Ubuntu 18.04` | `GCC-10`   | `Off`   | `Official`
| `test:x86-64:linux:gcc-10:cxx11-off:unsupported`    | `x86-64`  | `Ubuntu 18.04` | `GCC-10`   | `Off`   | `Unsupported`
| `test:x86-64:linux:gcc-10:cxx11-on:official`        | `x86-64`  | `Ubuntu 18.04` | `GCC-10`   | `On`    | `Official`
| `test:x86-64:linux:gcc-10:cxx11-on:unsupported`     | `x86-64`  | `Ubuntu 18.04` | `GCC-10`   | `On`    | `Unsupported`
| `test:x86-64:linux:clang-10:cxx11-off:official`     | `x86-64`  | `Ubuntu 18.04` | `Clang-10` | `Off`   | `Official`
| `test:x86-64:linux:clang-10:cxx11-off:unsupported`  | `x86-64`  | `Ubuntu 18.04` | `Clang-10` | `Off`   | `Unsupported`
| `test:x86-64:linux:clang-10:cxx11-on:official`      | `x86-64`  | `Ubuntu 18.04` | `Clang-10` | `On`    | `Official`
| `test:x86-64:linux:clang-10:cxx11-on:unsupported`   | `x86-64`  | `Ubuntu 18.04` | `Clang-10` | `On`    | `Unsupported`
| `test:aarch64:linux:gcc-10:cxx11-off:official`      | `AArch64` | `Ubuntu 18.04` | `GCC-10`   | `Off`   | `Official`
| `test:aarch64:linux:gcc-10:cxx11-off:unsupported`   | `AArch64` | `Ubuntu 18.04` | `GCC-10`   | `Off`   | `Unsupported`
| `test:aarch64:linux:gcc-10:cxx11-on:official`       | `AArch64` | `Ubuntu 18.04` | `GCC-10`   | `On`    | `Official`
| `test:aarch64:linux:gcc-10:cxx11-on:unsupported`    | `AArch64` | `Ubuntu 18.04` | `GCC-10`   | `On`    | `Unsupported`
| `test:aarch64:linux:clang-10:cxx11-off:official`    | `AArch64` | `Ubuntu 18.04` | `Clang-10` | `Off`   | `Official`
| `test:aarch64:linux:clang-10:cxx11-off:unsupported` | `AArch64` | `Ubuntu 18.04` | `Clang-10` | `Off`   | `Unsupported`
| `test:aarch64:linux:clang-10:cxx11-on:official`     | `AArch64` | `Ubuntu 18.04` | `Clang-10` | `On`    | `Official`
| `test:aarch64:linux:clang-10:cxx11-on:unsupported`  | `AArch64` | `Ubuntu 18.04` | `Clang-10` | `On`    | `Unsupported`
