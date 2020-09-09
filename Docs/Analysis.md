# astcenc 2.0 analysis

This page is a working document explaining my current understanding of the
performance and quality impact of some of the compression trials, refinement
passes, and the heuristics which drive them. To make things easier to
cross-reference, I've inserted some comments in the code to identify passes
that I've analyzed.

This information will hopefully allow as "tuned up" a future version of the
codec. At the moment it's just research.

## Configuration sensitivity

The most important thing to know for this analysis is that the impact of each
trial and refinement can vary dramatically depending on the codec
configuration. Principally this means:

* block size,
* quality preset,
* image color format.

For example, removing the initial 1 partition 1 plane fast-early-out
optimization actually improves the performance of 4x4 blocks in `-fast` and
`-medium`. However, it makes 5x5 blocks marginally slower and the slow down
increases as block size increases.

