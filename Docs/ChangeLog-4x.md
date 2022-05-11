# 4.x series change log

This page summarizes the major functional and performance changes in each
release of the 4.x series.

All performance data on this page is measured on an Intel Core i5-9600K
clocked at 4.2 GHz, running `astcenc` using AVX2 and 6 threads.

<!-- ---------------------------------------------------------------------- -->
## 4.0

**Status:** In development

The 4.0 release introduces some major performance enhancement, and a number
of larger changes to the heuristics used in the codec to find a more effective
cost:quality trade off.

* **General:**
  * **Feature:** The command line tool now has `-repeats <count>` for testing.
  * **Optimization:** Angular weight endpoint selection is restricted to
    QUANT_11 or lower. Higher quantization levels assume default 0-1 range,
    which is less accurate but faster.
  * **Optimization:** Maximum weight quantization for later trials is selected
    based on the best encoded weight quantization for the 1 plane 1 partition
    trial. This significantly reduces the search space for higher quant levels.

- - -

_Copyright © 2022, Arm Limited and contributors. All rights reserved._
