# symba
A market simulator built in C++.

## Building symba
To build Symba you will need [CMake](https://cmake.org/) 3.17+ and any C++ compiler with support for C++11. This project recommends [Ninja](https://ninja-build.org/) as the CMake build system.

This project requires libraries [GSL](https://www.gnu.org/software/gsl/) and [CLI11](https://github.com/CLIUtils/CLI11). They are available in both [vcpkg](https://github.com/microsoft/vcpkg) and [Conan](https://conan.io/) package managers.

<details><summary>Example: building on Windows with vcpkg</summary><p>

The following example builds Symba on Windows (cmd) using vcpkg and Ninja:

1. If your compiler is from Visual Studio, activate the x64 developer command prompt:
```
"<path-to-visualstudio>/VC/Auxiliary/Build/vcvars64.bat"
```

2. Install GSL using vcpkg:
```
"<path-to-vcpkg>/vcpkg" install gsl:x64-windows-static
```

3. Clone this project:
```
git clone https://github.com/andreasxp/symba
cd symba
```

4. Configure CMake in the build folder:
```
cmake -S . -B build -G "Ninja Multi-Config"
  -D CMAKE_CONFIGURATION_TYPES="Debug;Release"
  -D VCPKG_TARGET_TRIPLET=x64-windows-static
  -D CMAKE_TOOLCHAIN_FILE="<path-to-vcpkg>/scripts/buildsystems/vcpkg.cmake"
```

5. Build and install the project in Release configuration
```
cmake --build build --config Release
cmake --install build --config Release --prefix install
```

The executable can then be found in the `install/bin` directory.

</p></details>
<br/>

## Running symba
Symba has a command-line interface where you can specify the output directory. A basic example of running it would be:
```
./symba --output-dir <path-to-directory>
```

More information is available with `--help`.
