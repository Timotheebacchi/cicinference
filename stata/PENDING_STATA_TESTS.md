# Pending Stata Tests

Stata is not installed or not available on `PATH` in the current Codex
environment. The tests below must be run on a machine with an installed and
licensed Stata executable before claiming any Stata implementation passes.

## Commands

Using Stata/MP:

```bash
stata-mp -b do stata/tests/test_phase3_all.do
stata-mp -b do stata/tests/test_nosplit.do
stata-mp -b do stata/tests/test_split.do
stata-mp -b do stata/tests/test_kde.do
stata-mp -b do stata/tests/test_bootstrap.do
stata-mp -b do stata/tests/benchmark.do
```

Using Stata/SE:

```bash
stata-se -b do stata/tests/test_phase3_all.do
stata-se -b do stata/tests/test_nosplit.do
stata-se -b do stata/tests/test_split.do
stata-se -b do stata/tests/test_kde.do
stata-se -b do stata/tests/test_bootstrap.do
stata-se -b do stata/tests/benchmark.do
```

Using an executable named `stata`:

```bash
stata -b do stata/tests/test_phase3_all.do
stata -b do stata/tests/test_nosplit.do
stata -b do stata/tests/test_split.do
stata -b do stata/tests/test_kde.do
stata -b do stata/tests/test_bootstrap.do
stata -b do stata/tests/benchmark.do
```

## Current Status

- R reference generation has been executed with `Rscript`.
- Static Stata/Mata inspection has been performed during migration work.
- No-split, split, and KDE reference tests exercise the public
  `cicinference` command.
- Bootstrap reference tests exercise the deterministic common-index bootstrap
  internals; live Stata bootstrap RNG output is not expected to equal R by seed.
- Actual Stata execution remains pending.
