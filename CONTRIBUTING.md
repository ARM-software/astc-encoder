# Contributing

Contributions are welcome.

## Licensing and copyright

By submitting a contribution to this project, you agree to license it under the
[Apache License 2.0](LICENSE.txt).

We do not require copyright assignment. You or your employer, as applicable,
retain the copyright in your contribution. If you wish, you may add your own
copyright notice to the header of any file you modify. Do not modify any other
part of the license text in the header.

## General policies

For minor changes and bug fixes, please raise a pull request directly. An issue
is not required.

For major changes, please raise an issue and discuss it with the project team
before starting work on a pull request. This lets us align on the idea and
approach before anyone spends significant time on changes that cannot be
accepted or that may need reworking later.

Each issue or pull request should cover only one feature, optimization, or bug
fix. This makes it much easier to manage discussion threads and validate
changes.

Pull requests should describe the motivation for the change, the tests
performed, and the architectures affected. Performance optimizations should
also include benchmark results that demonstrate the impact of the change.
Include enough details so that reviewers can reproduce your results.

Pull requests must follow the style of the surrounding code, and ensure that
builds remain free of compiler warnings.

## Raising bug reports

When reporting a bug, please include:

- The astcenc version or git revision used.
- The compiler and compiler version used.
- The CMake build configuration used.
- The complete command line used, if using the command line tool.
- A reproducer image, if you are able to share it.

## Codec image quality

The codec uses strict floating-point behavior and is designed to produce
invariant results across compilers and platforms. The test suite checks that
pull requests produce bit-exact results. Small image-quality differences
introduced by a pull request, such as those resulting from arithmetic
reassociation, are usually acceptable, but require the test reference results
to be rebuilt.

Our default policy is that changes should not cause image-quality regressions
beyond differences caused by arithmetic reassociation. Major changes that shift
the performance-quality balance should be discussed in an issue before starting
work on a pull request.

## Codec performance

Performance optimizations must not regress performance on any supported
architecture. Changes that improve some architectures while leaving others
unchanged are acceptable.

## Codec use of SIMD intrinsics

We encapsulate architecture-specific behavior behind the
`astcenc_vecmathlib_<arch>.h` vector library abstraction. Code outside this
abstraction must be architecture-agnostic.

If an ISA provides a useful operation that requires extending the vector
library interface, the operation must be implemented for all supported
architectures, using multi-instruction emulation where necessary. This
restriction reduces some optimization opportunities but makes the code much
easier to maintain.

- - -

_Copyright © 2026, Arm Limited and contributors._
