empi
====

University of Warsaw, Department of Biomedical Physics ⓒ 2015–2016  
Enhanced Matching Pursuit Implementation (empi)  
Author: Piotr Różański <piotr@develancer.pl>  
& improvements of code and build process thanks to Aleks Chrabrow

## What is empi?

empi is an implementation of Matching Pursuit algorithm
([Mallat, Zhang 1993](http://dx.doi.org/10.1109/78.258082))
with optimal Gabor dictionaries
([Kuś, Różański, Durka 2013](http://dx.doi.org/10.1186/1475-925X-12-94)).
It is a highly-optimized multithreaded version written in C++,
designed as a faster replacement for
[MP5](https://github.com/BrainTech/matching-pursuit). Therefore, it shares most
of the input/output specification with MP5, and can be used as MP decomposition
tool in [SVAROG](https://github.com/BrainTech/svarog).

The goal is an optimal decomposition of the input signal as a linear combination
of functions from predefined set (dictionary) of Gabor atoms, including ordinary
Gaussians as a special case for frequency = 0.

## How to get empi?

You can compile empi from source, or download the precompiled versions from the
“Releases” tab. Both are available on
[project's GitHub](https://github.com/develancer/empi). If you decide to use
the precompiled binaries, you can skip the “Compilation” section altogether.

### Compilation

#### Requirements

To compile empi, CMake build system is required. The only external library
requirement is the FFTW library in version 3. Both library
and the development headers must be installed to compile empi. Under Ubuntu,
package “libfftw3-dev” does the trick. MacOS and other Linux distributions
may have packages of slightly different names.

Also, you will need a modern C++ compiler with support for C++11 standard.
OpenMP support is recommended.

#### Configuration

This project uses CMake, so the proper way of compilation depends on your
environment. Generally speaking, you can use CMake-gui to generate build files
for your specific configuration.

#### Unix-family OS

The easiest way is to run

	cmake .

_or_, if you need to build standalone binaries
(it will disable some platform-specific optimizations),

	cmake -DSTANDALONE=1 .

followed by

	make

in the directory where you cloned your repository; you can also do an
out-of-source build, if you prefer. If successful, binary file “empi” shall
appear. It can by copied to /usr/local/bin (in order to be available in $PATH)
by calling

	sudo make install

#### Cross-compilation

Binary versions for 32-bit or 64-bit Microsoft Windows can be cross-compiled
under Linux. Under Ubuntu, packages _g++-mingw-w64-i686_ (for 32-bit) and
_g++-mingw-w64-x86-64_ (for 64-bit) are needed. In different distributions,
package names may vary. Also, shared versions of FFTW
libraries (DLL and include files) should be already installed in each toolchain.
To generate cross-compiled exe files, execute

	cmake -DSTANDALONE=1 .
	make empi-win32.exe empi-win64.exe

To use generated binaries under MS Windows, the dynamic version of FFTW library
(DLL file) has to be placed in the same directory as the executable file.
All the other libraries will be linked statically into the executable.

## How to use empi?

Single invocation of empi will

* read a single binary file (or its part),
* decompose it as a linear combination of Gabor atoms, and
* save the results to a specified format (currently only
[SVAROG](https://github.com/BrainTech/svarog)'s “book” format is supported).

empi needs to be run with a single command-line argument: a path to the
configuration file. If run with no arguments, it will print the correct usage.
For backward compatibility with MP5, all arguments starting with “`-`” are ignored.

### Configuration file format

Let us start with a sample configuration file:

	energyError 0.01
	maximalNumberOfIterations 50
	energyPercent 99.0

	MP SMP

	nameOfDataFile signal.bin
	nameOfOutputDirectory .
	samplingFrequency 128.0

	numberOfChannels 3
	selectedChannels 1-3

Most parameters are straightforward, but we shall describe them one by one:

* _energyError_ is the ε² parameter in optimal Gabor dictionary construction.
Usually the values will be close to 0. Smaller value will allow for a more
precise decomposition, but it will also engage more time and RAM.

* _maximalNumberOfIterations_ is the upper limit for the number of iterations,
and therefore, a maximal number of atoms in the resulting decomposition.

* _energyPercent_ is the percent of the energy that we require to be “explained”
by decomposition, Requiring 99% of energy to be explained means that we will be
performing decomposition until residual energy fall below 1% of the total energy
of the signal.

The decomposition will iterate until _maximalNumberOfIterations_ or
_energyPercent_ will be fulfilled, whichever comes first.

* _MP_ is a selected variant of Matching Pursuit. Currently, only SMP is supported.

* _nameOfDataFile_ is a path to the binary signal file, relative to the current
directory. The input file should consist of 32-bit float values in the byte order
of the current machine (no byte-order conversion is performed). For multichannel
signals, first come the samples for all channels at t=0, then for all channels
at t=Δt, and so forth. In other words, the signal should be written
in column-major order (rows = channels, columns = samples).

* _nameOfOutputDirectory_ is a path to the output directory, relative to the
current directory. The output file will be named based on the name of the input
file, e.g. if input file is “signal.bin”, the output will be named “signal_XYZ.b”,
where XYZ is the selected variant of MP (parametr _MP_).

* _samplingFrequency_ is a sampling frequency of the input signal, specified in
hertz.

* _numberOfChannels_ is a number of all channels in the input signal.

* _selectedChannels_ specify which channels should be read from the signal and
decomposed. These can be specified as a single channel `1`, as an interval `1-5`,
as a list `3,4` or mixed: `1-2,5,8-10`.

## Disclaimer

empi is free software; you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

empi is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
empi (file “LICENCE”); if not, write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
