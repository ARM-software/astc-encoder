# Branch 2.0 change log

This page maintains an overview of the major changes made on the 2.0 branch of
the compressor, both in terms of algorithm and implementation. To help make
bisection easier, the list is maintained in approximate chronological order
(oldest first).

## Maths library changes

The astc_mathlib library has been cleaned up, with all unused functionality
stripped out.

The vector template classes have been stripped back to a basic implementation,
with the flexible channel swizzle functionality removed. This makes it easier
to see what is going on in the code (less magic), and dramatically improves
compile time due to the reduction in template processing, at the expense of
a little more scalar code expansion where swizzles could have been used.

## Implementation changes

