# Branch 2.0 design notes

# Objectives

The aim of the 2.0 work is to achieve the equivalent of `-thorough` encoding
quality with `-medium` performance. This equates to an ~8x performance
improvement over the version 1.x codec, at the same output quality, and without
losing any of the format features (e.g. block sizes and modes) currently
supported.

There are four main approaches which will be used to achieve this.

* Code simplification
* Datapath restructuring
* Algorithm tuning
* CPU SIMD vectorization

We are currently focussing on the code simplification and algorithm review.

## Code simplification

The version 1.0 codec was initially developed as part of the the ASTC format
design, evolving as the format was developed and refined. A lot of the v1.x
implementation was therefore built to be modular and easy to adapt as the
format design changed, rather than for the best possible performance. Now that
the format has been published we can review all of the existing code, and
rework it for performance rather than ease of modification.

## Datapath restructuring

The current code does one block at a time, in a depth-first search with full
iterative refinement at each stage until the target quality is either hit or we
run out of codec options to try. This approach is easier to understand, but the
resulting code is control plane heavy and includes a lot of overhead which
cannot be amoritzed over multiple blocks. As a first step for SIMD
vectorization the data path must be reviewed to make sure that it can be
vectorized efficiently.

* Batching processing passes for multiple blocks.
* Reducing branching and branch divergence to allow SIMD batching.
* Reducing scatter/gather loads to ensure efficient memory access.
* Amortizing memory access cost over multiple blocks

The restructure becomes particular important with an eye on a long-term
objective (v3?) to support compression passes running on a GPU, as it is
necessary to offload relatively large batches of work, while also avoiding
bottlenecking GPU performance on memory access overheads rather than useful
processing.

## Algorithm tuning

As part of the restructuring we will review the effectiveness of the existing
optimization passes that we have in the code, to see what can be removed. We
know that some passes are expensive and/or difficult to vectorize, and do not
actually justify that cost in significant quality improvements.

In addition many math library operations will be reviewed and replaced with
faster implementations which are less accurate but more friendly towards
vectorization in future.

## CPU SIMD vectorization

The final stage, once we are happy with the data path design, will be to look
at implementing a true SIMD data path for the CPU. This could be done using
manual SIMD code, such as intrinsics, or a compiled SIMD generator,
such as ISPC; the decision here can happen later.

# TODOs

* The `double4` vector type doesn't vectorize well with SIMD; it should be
  removed, but this probably means we need a new algorithm for that part
  of the code.
