# Performance Notes

Status: benchmark script written, actual Stata timing pending.

## Expected Complexity

- `nosplit`: sorting plus binary-search rectangle counts, roughly `O(n log n)`.
- `split`: same order as `nosplit`, using deterministic first/last half-samples.
- `kde`: sorted centered `Y`, prefix sums, and binary searches, roughly `O(n log n)`.
- `bse` and `bpc`: `O(B * n log n)` in the current readable Mata implementation.

## Pending Benchmark Command

Run from the repository root on a machine with Stata:

```bash
stata-mp -b do stata/tests/benchmark.do
```

Equivalent local executable names may be `stata-se` or `stata`.

The current Codex environment does not have Stata on `PATH`, so no Stata
runtime measurements are claimed yet.
