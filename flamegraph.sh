#!/bin/bash -ex
perf script | ../FlameGraph/stackcollapse-perf.pl > out.perf-folded
../FlameGraph/flamegraph.pl --title=flamegraph out.perf-folded > flamegraph.svg
