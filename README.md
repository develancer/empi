empi
====

Enhanced Matching Pursuit Implementation (empi)  
Author: Piotr Różański <piotr@develancer.pl> ⓒ 2015–2020

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
However, since the purpose of the provided binaries is to be as compatible as
possible, they may not take full advantage of your specific architecture. To
achieve maximal performance, compiling empi from source is recommended.

Additionally, precompiled binaries for OS X does not support OpenMP (therefore,
they run a single thread). To use OpenMP on OS X, compilation is necessary.

### Compilation

#### Requirements

To compile empi, CMake build system is required. The only external library
requirement is the FFTW library in version 3. Both library
and the development headers must be installed to compile empi. Under Ubuntu,
package “libfftw3-dev” does the trick. MacOS and other Linux distributions
may have packages of slightly different names. Under Windows, follow the
FFTW installation instructions.

Also, you will need a modern C++ compiler with support for C++11 standard.
OpenMP support is recommended. Under OS X, OpenMP runtime
(available at [LLVM Download Page](http://llvm.org/releases/download.html))
should be downloaded
and put into system library directory (e.g. /usr/lib) prior to compiling.

#### Configuration

This project uses CMake, so the proper way of compilation depends on your
environment. Generally speaking, you can use CMake-gui to generate build files
for your specific configuration.

#### Unix-family OS

The easiest way is to run

	cmake .

_or_, if you need to build standalone binaries
(however, it will disable some platform-specific optimizations),

	cmake -DSTANDALONE=1 .

followed by

	make

in the directory where you cloned your repository; you can also do an
out-of-source build, if you prefer. If successful, binary file “empi” shall
appear. It can by installed to system directory (e.g. /usr/local/bin)
by calling

	sudo make install

#### Cross-compilation

Binaries can also be cross-compiled for a different system or architecture.
Detailed information on cross-compiling empi can be found in the Appendix,
at the end of this document.

## How to use empi?

Single invocation of empi will

* read a single binary file (or its part),
* decompose it as a linear combination of Gabor atoms, and
* save the results to either SQLite, JSON or legacy format
([SVAROG](https://github.com/BrainTech/svarog)'s “book” format is supported).

empi needs to be run with a single command-line argument: a path to the
configuration file. If run with no arguments, it will print the correct usage.
Default output format is SQLite; however:
* if `-x` is given as an argument, legacy “book” file will be created instead;
* if `-j` is given as an argument, JSON output file will be created instead.

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

	numberOfSamplesInEpoch 1280
	selectedEpochs 1-2

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

* _MP_ is a selected variant of Matching Pursuit. Following variants are supported:

	* “SMP” decomposes every channel independently

	* “MMP1” finds, in every iteration, set of atoms which differ only in amplitude,
	optimizing sum of squares of the scalar products across channels

	* “MMP2” finds, in every iteration, set of atoms which differ only in amplitude,
	optimizing sum of the scalar products across channels (much faster than MMP1)

	* “MMP3” finds, in every iteration, set of atoms which differ in amplitude
	and/or phase, optimizing sum of squares of the scalar products across channels

* _nameOfDataFile_ is a path to the binary signal file, relative to the current
directory. The input file should consist of 32-bit float values in the byte order
of the current machine (no byte-order conversion is performed). For multichannel
signals, first come the samples for all channels at t=0, then for all channels
at t=Δt, and so forth. In other words, the signal should be written
in column-major order (rows = channels, columns = samples).

* _nameOfOutputDirectory_ is a path to the output directory, relative to the
current directory. The output file will be named based on the name of the input
file, e.g. if input file is “signal.bin”, the output will be named either
“signal_XYZ.b” or “signal_XYZ.json” (depending on the selected output format),
where XYZ is the selected variant of MP (parametr _MP_).

* _samplingFrequency_ is a sampling frequency of the input signal, specified in
hertz.

* _numberOfChannels_ is a number of all channels in the input signal.

* _selectedChannels_ specify which channels should be read from the signal and
decomposed. These can be specified as a single channel `1`, as an interval `1-5`,
as a list `3,4` or mixed: `1-2,5,8-10`.

* _numberOfSamplesInEpoch_ specifies the size of each signal segment (in samples);
segments will be processed in order, and their decomposition will be written
to the same output file. This parameter is optional; if not given, the entire
signal will be processed as a single segment.

* _selectedEpochs_ is only valid with _numberOfSamplesInEpoch_ and it specifies
a list of epoch numbers (starting from 1) to be processed. These can be passed
as a single epoch `1`, as an interval `1-100`, as a list `1,2,3` or mixed:
`1-100,201-300,400`. This parameter is optional; if not given, then all epochs
(the entire signal) will be processed.

* _minAtomScale_ specifies the minimum scale of atoms in the MP dictionary.
If not specified, it defaults to the scale length of 2 samples
(scale in seconds = 2 / sampling frequency)

* _maxAtomScale_ specifies the maximum scale of atoms in the MP dictionary.
If not specified, it defaults to the length of the analysed signal.

* _maxAtomFrequency_ specifies the maximum frequency (in hertz)
of atoms in the MP dictionary.

Specifying _minAtomScale_, _maxAtomScale_ and _maxAtomFrequency_ is not
mandatory, but it allows to greatly reduce the computational time.

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

## Appendix. Cross-compiling empi

This section is dedicated to the steps needed to cross-compile empi for a
different operating system and/or architecture on 64-bit Linux
(specifically Ubuntu 15.10).

### Apple Mac (OS X)

This is probably the most complex case. The project
[osxcross](https://github.com/tpoechtrager/osxcross) is a good start. It allows
to build a full cross-compiling toolchain for both 32-bit and 64-bit OS X
compilation. This section will focus on 64-bit toolchain (32-bit is analogous).
As a part of its configuration, it is necessary to sign up for the Apple
developer's account and download the Xcode package from the official Apple site
(don't worry, it's free).

After the toolchain is built, it is necessary to build the FFTW library
and copy it (the library itself and header files) into the toolchain.
To successfully compile FFTW, it is necessary to specify all the paths to
the compiler binaries, e.g. (some of the below may not be necessary)

	AR=x86_64-apple-darwin15-ar \
	AS=x86_64-apple-darwin15-as \
	CC=x86_64-apple-darwin15-clang \
	LD=x86_64-apple-darwin15-ld \
	RANLIB=x86_64-apple-darwin15-ranlib \
	./configure --host=x86_64-apple
	make

The static version of compiled FFTW library (libfftw3.a) has to be placed in
the toolchain library directory.

Having it all set, return to the empi directory, run

	make empi-osx64

and that's it, finally. The cross-compiled version will _not_ use OpenMP.

### Microsoft Windows

Both 32-bit or 64-bit Windows binaries can be cross-compiled. Packages
_g++-mingw-w64-i686_ (for 32-bit) and _g++-mingw-w64-x86-64_ (for 64-bit)
must be installed. Since FFTW developers generally discourage manually
compiling FFTW for Windows, it is better to stick with shared versions,
which can be obtained from
[FFTW download page](http://www.fftw.org/install/windows.html).
After downloading, install DLL files and include headers in each toolchain
(i686 and/or x86-64), e.g.:

* /usr/lib/gcc/i686-w64-mingw32/5.3-win32 for 32-bit Windows
* /usr/lib/gcc/x86_64-w64-mingw32/5.3-win32 for 64-bit Windows

To generate cross-compiled exe files, execute

	make empi-win32.exe empi-win64.exe

To use generated binaries under MS Windows, the dynamic version of FFTW library
(DLL file) has to be placed in the same directory as the executable file.
All the other libraries will be linked statically into the executable.
