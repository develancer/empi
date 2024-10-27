empi
====

Enhanced Matching Pursuit Implementation (empi)  
Author: Piotr Różański <piotr@develancer.pl> ⓒ 2015–2023

**Cite as:**

> Piotr T. Różański:
> _empi_: GPU-Accelerated Matching Pursuit with Continuous Dictionaries,  
> ACM Transactions on Mathematical Software, Volume 50, Issue 3, Article 17 (2024),
> DOI: [10.1145/3674832](https://doi.org/10.1145/3674832)

## What is empi?

empi is an implementation of Matching Pursuit algorithm
([Mallat, Zhang 1993](http://dx.doi.org/10.1109/78.258082))
with optimal dictionaries
([Kuś, Różański, Durka 2013](http://doi.org/10.1186/1475-925X-12-94))
including simulation of continuous (quasi-infinite) dictionaries
([Różański 2024](https://doi.org/10.1145/3674832)),
supporting both Gabor atoms as well as atoms related to non-Gaussian envelopes
([Różański 2020](https://doi.org/10.1049/iet-spr.2019.0246)).
It is a highly-optimized multi-threaded version written in C++, with GPU support,
designed as a faster replacement for
[MP5](https://github.com/BrainTech/matching-pursuit).
Additionally, this implementation can be used as MP decomposition
tool in [SVAROG](https://github.com/BrainTech/svarog)’s most
[up-to-date version](https://gitlab.com/fuw_software/svarog2-packager/-/releases).

The goal is to provide an optimal decomposition of the input signal as a linear
combination of functions from predefined set (dictionary) consisting mainly of
oscillating atoms. By combining the optimal dictionary construction with a
detailed analysis of the maximum error within a single iteration,
it can be used to simulate a continuous dictionary,
therefore relieving the user from the necessity of defining the particular
structure for the dictionary. What is even more important, it therefore completely
eliminates the statistical bias caused by the dictionary structure—an important
effect which has been earlier dealt with by e.g. introduction of stochastic
dictionaries.

There are two modes of CPU parallelization, which can also be used together.
First, a number of independent workers can be started, and each worker will
process a separate subset of segments and/or channels. Second, each worker can
be started with a number of concurrent CPU threads. Either way, in addition to
all workers’ threads, an additional single thread will be active and responsible
for writing the decomposition results from all workers, in proper order, to the
output file.

If the GPU devices are used in addition to CPU, each worker is started with one
additional thread (corresponding to a separate CUDA stream) for each GPU device.
The performance gain from enabling GPU devices, especially those with good
double-precision floating point capabilities, is significant.

_empi_ includes code from CLI11 and SQLite projects in the “vendor” subdirectory.

## How to get empi?

You can compile empi from source, or download the precompiled versions from the
“Releases” tab. Both are available on
[project’s GitHub](https://github.com/develancer/empi). If you decide to use
the precompiled binaries, you can skip the “Compilation” section altogether.
However, since the purpose of the provided binaries is to be as compatible as
possible, they may not take full advantage of your specific architecture. To
achieve maximal performance and/or use GPU in calculations, compiling empi
from source is recommended.

### Compilation

#### Requirements

To compile empi, CMake build system is required. The only external library
requirement is the FFTW library in version 3. Both library
and the development headers must be installed to compile empi. Under Ubuntu,
package “libfftw3-dev” does the trick. MacOS and other Linux distributions
may have packages of slightly different names. Under Windows, follow the
FFTW installation instructions.

Also, you will need a modern C++ compiler with support for C++17 standard,
as well as CMake version 3.12 or later.

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

	make empi

in the directory where you cloned your repository; you can also do an
out-of-source build, if you prefer. If successful, binary file “empi” shall
appear. It can be installed to system directory (e.g. /usr/local/bin)
by calling

	sudo make install

## How to use empi?

Single invocation of empi will

* read a single binary signal file (or its part),
* decompose it as a linear combination of well-defined structures, and
* save the results as either SQLite (default) or JSON

Directory _demo_ includes Python and Matlab/Octave scripts demonstrating
how to access data from the resulting SQLite decomposition file.

empi needs to be run with at least two command-line argument: a path to the input
file and a path to the output file. It can be run with `--help` flag to list all
possible flags and arguments:

```
Enhanced Matching Pursuit Implementation (empi) 1.0.0
Usage: empi [OPTIONS] input_file output_file

Positionals:
  input_file TEXT REQUIRED    Path to the input signal file or input configuration file
  output_file TEXT REQUIRED   Path for the output file unless configuration file is used

Options:
  -h,--help                   Print this help message and exit
  --version                   Print version and exit
  -c INT=1                    Number of channels in the input signal
  -f FLOAT                    Sampling frequency of the input signal in hertz (default: 1 Hz)
  -i INT                      Maximum number of iterations (default: no limit)
  -o TEXT                     Parameter optimization mode: none|local|global (default: global)
  -r FLOAT=0.01               Energy of the residual as a fraction of the total signal energy
  --channels TEXT             Range of channels to process, e.g. 1-3,5,8-9 (default: all)
  --cpu-threads UINT=6        Number of CPU threads for each worker
  --cpu-workers UINT=1        Number of independent CPU workers to run
  --delta                     Include delta-type atoms
  --full-atoms-in-signal      Prohibit atoms from exceeding the time range of the signal
  --energy-error FLOAT=0.05   Epsilon-squared parameter corresponding to the dictionary size
  --dictionary-output TEXT    Path to create a dictionary structure XML file (default: none)
  --gpu-id TEXT               Comma-separated ID list of GPU device(s) to use (default: none)
  --input64                   Read input data as double-precision (64-bit) floating point values (default: read as 32-bit values)
  --mmp1 Excludes: --mmp3     Use multi-variate decomposition with constant phase across channels
  --mmp3 Excludes: --mmp1     Use multi-variate decomposition with variable phase across channels
  --opt-max-iter INT=10000    Maximum number of iterations for local parameter optimization
  --opt-target FLOAT=1e-05    Target accuracy (relative to the initial dictionary size) for local parameter optimization
  --residual-log-dir TEXT     Directory in which residual energy log files should be created (default: none)
  --segment-size INT          Number of samples in each segment (default: all samples)
  --segments TEXT Needs: --segment-size
                              Range of signal segments, e.g. 1-100,201-300 (default: all)
  --gabor                     Include atoms with Gaussian envelope (not needed if any other --gabor-* option is given)
  --gabor-freq-max FLOAT      Maximum frequency (in hertz) for Gaussian envelope (default: auto)
  --gabor-scale-min FLOAT     Minimum scale (in seconds) for Gaussian envelope (default: auto)
  --gabor-scale-max FLOAT     Maximum scale (in seconds) for Gaussian envelope (default: auto)
  --gabor-half-width FLOAT=1.5
                              Half-width of the Gaussian envelope function
```

### Command-line options

There are two required positional arguments:
1. _input_file_ is the full (or relative to the current directory) path to the input file.
The input file should consist of 32-bit (or 64-bit if the `--input64` flag is given)
floating-point values in the byte order  of the current machine
(no byte-order conversion is performed). For multichannel
signals, first come the samples for all channels at t=0, then for all channels
at t=Δt, and so forth. In other words, the signal should be written
in column-major order (rows = channels, columns = samples).
2. _output_file_ is the full (or relative to the current directory) path for the output file.
If the path ends in `.json`, JSON-formatted text file will be created. 
Otherwise, SQLite database file will be created.

The optional parameters are described below:

#### Properties of the input file

* `-c` represents the number of all channels in the input file, the default corresponding to a single-channel signal (as with `-c 1`).
* `-f` represents the signal's sampling frequency in hertz, the default being 1 Hz.
* `--channels` allows to specify the subset of channels (between 1 and the value of `-c`)
that should be read from the signal and decomposed.
These can be specified as a single channel `1`, as an interval `1-5`, as a list `3,4` or mixed: `1-2,5,8-10`.
* `--input64` assumes the input signal file consists of 64-bit floating point values, as opposed to the 32-bit as default.
* `--segment-size` specifies the size of each signal segment (in samples);
  segments will be processed independently, and their decomposition will be written
  to the same output file. If this parameter is absent, the entire
  signal will be processed as a single segment.
* `--segments` is only valid with `--segment-size`, and it specifies
a list of epoch numbers (starting from 1) to be processed. These can be passed
as a single epoch `1`, as an interval `1-100`, as a list `1,2,3` or mixed:
`1-100,201-300,400`. If not given, all epochs (the entire signal) will be processed.

#### Decomposition process

The decomposition will iterate until `-i` iterations or
`-r` residual energy is reached, whichever comes first.

* `-i` specifies an upper limit for the number of iterations,
  and therefore, a maximal number of atoms in the resulting decomposition.
* `-r` is the percent of the residual energy that can be left un-explained
  by decomposition. For example, specifying `-r 0.01` corresponds to
  performing the decomposition until the energy of the residual
  falls below 1% of the initial energy of the signal.
* `--dictionary-output` value allows to generate an MPTK-style XML file with dictionary structure at a given path.
* `--mmp1` specifies a constant-phase multi-variate decomposition,
  while `--mmp3` specifies a variable-phase variant. If neither is given, 
  each channel is processed separately. Detailed description of multi-variate
  decomposition modes can be found in the literature, e.g.
  [Kuś, Różański, Durka 2013](http://doi.org/10.1186/1475-925X-12-94).

#### Structure of the dictionary

* `-o` specifies how the local optimization will be used in a parameter space.
`-o global` is default and results in a full simulation of a continuous dictionary.
`-o local` enables local parameter optimization, but starting only from each iteration’s best match
and therefore, does not guarantee choosing the globally best atom in each iteration. However, it could constitute
a sensible tradeoff between precision and performance for real-world usages, as it can be much faster than `-o global`.
Specifying `-o none` results in using only a discrete dictionary which structure will be specified by
the `--energy-error` flag as described below. In this case, one should also request an `--energy-error` value
smaller than the default of 0.05 (e.g. 0.01).

* `--delta` enables the use of delta atoms in the dictionary.

* `--energy-error` specifies the ε² parameter in optimal dictionary construction.
Usually the values will be close to 0. Smaller value will allow for a more
precise decomposition, but it will also engage more time and RAM.
When simulating the continuous dictionary (`-o` option), this parameter does not affect
the results, only time and memory consumption. The default value of 0.05 is approximately
optimal for simulating continuous dictionaries, but for `-o none` one should request a smaller value (e.g. 0.01).

* `--full-atoms-in-signal` restricts the possible positions and scales of atoms in the dictionary
so that they fully overlap with the signal. This option should be set to obtain full compatibility with MPTK.
Otherwise (default), the atoms may exceed the signal range
and the signal will be assumed to have zero values outside the actual range.

* `--gabor` enables the use of Gabor atoms in the dictionary.
This flag is implicit if any of the `--gabor-*` parameters is set.

* `--gabor-scale-min` specifies the minimum scale of Gabor atoms in the dictionary.
If not specified, the minimum scale is taken as the shortest scale permitted
by the dictionary construction for the given value of ε².

* `--gabor-scale-max` specifies the maximum scale of Gabor atoms in the dictionary.
If not specified, it defaults to the length of the signal segment.

* `--gabor-freq-max` specifies the maximum frequency (in hertz)
of Gabor atoms in the dictionary. If not specified, defaults to Nyquist frequency.

Specifying `--gabor-scale-min`, `--gabor-scale-max` and `--gabor-freq-max` is not
mandatory, but it allows to additionally reduce the computational time.

#### Parallelization scheme

* `--cpu-workers` defines the number of independent workers, where each
worker will independently analyse different subset of signal segments.
* `--cpu-threads` defines the number of CPU computation threads per each worker.
* `--gpu-id` (only if compiled with GPU support) defines the comma-separated list
of GPU devices that should assist in the decomposition.

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

### Disclaimer for the included source code of the CLI11 project

CLI11 2.1.2 Copyright (c) 2017-2021 University of Cincinnati, developed by Henry
Schreiner under NSF AWARD 1414736. All rights reserved.

Redistribution and use in source and binary forms of CLI11, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
