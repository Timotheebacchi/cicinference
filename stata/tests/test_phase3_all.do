version 16.0
clear all
set more off

do "stata/tests/test_phase3_rank_quantile.do"
do "stata/tests/test_phase3_nosplit_deterministic.do"
do "stata/tests/test_phase3b_density_variance.do"

display as result "All Phase 3 utility tests passed."
