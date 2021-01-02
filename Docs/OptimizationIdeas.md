# Optimization ideas

This page contains a relatively raw dump of optimization ideas we've had, and
areas we want to investigate.

Most of these are raw ideas that have not been entirely thought thorough. Some
will no doubt turn out to be rubbish ideas that don't work, but hopefully there
are a few gems in here =)

## Optimizations

This section lists specific optimization ideas.

### Postpone full iterative refinement

Iterative refinement is expensive, and not always beneficial. We currently
refine every trial based on estimated error, and then keep the best based on
actual error.

To reduce the amount of refinement done, but without losing too much quality,
we should see if it is possible to reduce the amount of refinement done before
candidate selection and then only refine the final one or two candidates (but
potentially refine them more than we do today).

### Heuristic based on block / weight decimation anisotropy

Decimated weight grids can end up quite anisotropic, especially at larger
block sizes. For example, the most extreme is a 2x12 grid for a 12x12 block.

Gut feel says we should be able to exploit this - either a straight-heuristic
to reject very unbalanced block sizes, or an intelligent heuristic which
looks at the variability in the input data in each axis to determine which
weight decimations to consider. Could also reject very high/low frequency
grids, if they are a bad fit, not just anisotropic ones.

This will mostly benefit the larger block sizes on `-medium` and `-thorough`
compression modes; percentiles tend to filter out the most unbalanced modes.

### Block mode and Decimation table sorting

The current code in `astcenc_block_sizes2.cpp` builds Block Modes based on the
value of the encoded mode index. The decoded block mode behavior is not nicely
sequential with either the likely utility of the block, or the use of data
resources. This gives poor locality and makes other optimizations harder.

- **Option one:** Sort block modes by the decimation entry that they reference,
  this will allow better spatial locality when processing sequential block
  modes that use the same table. These could be sorted by weight count so that
  there is some form of "ladder" in terms of weight count bitrate usage.
- **Option two:** Sort block modes by their "usefulness" centile, and test the
  most valuable ones first. This will allow us to build an early-out
  mechanism which tests the most useful block modes first. (Today we cut
  statically, based on the quality preset, but there is no dynamic trimming).
