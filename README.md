# RiskOddsCalculator

Calculates the odds of winning a battle in the game Risk
based on Markov chain probabilities, as described in this paper:
http://www4.stat.ncsu.edu/~jaosborn/research/RISK.pdf

Requires jblas library for linear algebra:
http://jblas.org
(written using V1.2.4)

Written by John Torgerson, 2016
Feel free to use or modify for any purpose.

The matrix operations seem to be too intensive to work with any
kind of large numbers, at least as written. It may be possible
to optimize things quite a bit, but I haven't tried. A lot of
this code was quick and dirty, just to get things working well
enough for my purposes, which this does.
