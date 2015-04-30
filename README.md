## Installation ##


```
#!r

library(devtools)
install_bitbucket("mkuhn/parallelRandomForest", ref="parallelRandomForest")

```

## Introduction ##
This repository contains a modified version of the original RF package for R.
The function randomForest behaves identical, with two additional parameters:
skip.checks (default: true) disables checking the matrices for NAs, and
nthreads specifies the number of threads to use.

In addition, "x" has to be a "raw" matrix, i.e. it can contain only single
bytes. As the default conversion can lead to unexpected results (for NAs and
values >255), it has to be done by the user.

The algorithm is optimized for "x" to contain only few unique values that
should be as small as possible. Ideally, all values in "x" are less or equal
than 2: in this case, a column only needs to be scanned twice to decide the
best split.

Using more than one thread creates multiple forests that are then combined,
due to the implementation of  the "combine" function this means that some
fields of the result are empty. Note that if your whole forest could be grown
in seconds, it doesn't make sense to use threads as there is some overhead (in
terms of time) to start a new subprocess.