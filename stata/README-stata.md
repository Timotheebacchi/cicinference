# cicinference Stata Migration

This directory contains the Stata port accompanying the R package
`quantcdf.inference`.

The GitHub source tree intentionally includes development material such as
reference CSVs, Stata `.do` tests, migration notes, and validation scripts.
Release archives for normal Stata users are built separately and exclude those
development files.

## Development Files

Tracked development files include:

- `stata/tests/`
- `stata/tests/reference_outputs/`
- `stata/MIGRATION_MAP.md`
- `stata/PENDING_STATA_TESTS.md`
- `stata/internal/`

These files are for migration, regression testing, reference generation, and
maintenance. Do not add them to `.gitignore`.

## Regenerate R References

From the repository root:

```bash
Rscript stata/tests/make_reference_outputs.R
```

This regenerates CSV inputs, R reference outputs, and intermediate utility
references readable from Stata. Bootstrap validation references include
R-exported common resampling indices because R and Stata RNGs differ.

## Run Stata Tests

These commands require an installed and licensed Stata executable:

```bash
stata-mp -b do stata/tests/test_phase3_all.do
stata-mp -b do stata/tests/test_nosplit.do
stata-mp -b do stata/tests/test_split.do
stata-mp -b do stata/tests/test_kde.do
stata-mp -b do stata/tests/test_bootstrap.do
stata-mp -b do stata/tests/benchmark.do
```

Equivalent executable names may be `stata` or `stata-se` depending on the local
installation.

The current development environment used by Codex does not have Stata on
`PATH`, so Stata tests are written but not claimed to pass until run locally.

## Build A User Release

From the repository root:

```bash
stata/build_release.sh
```

The script creates `stata/release/`, copies only user-facing runtime and
documentation files, and creates `cicinference-stata-release.zip`.

The generated `stata/release/` directory is ignored by Git. Development files
such as tests, CSV references, migration notes, logs, and benchmarks are not
copied into the release directory.

## Files Normal Users Receive

A release archive is intended to contain only:

- `cicinference.ado`
- `cicinference.sthlp`
- runtime Mata source files required by `cicinference.ado`
- `cicinference.pkg`
- `stata.toc`

GitHub source clones contain the full development test suite and R reference
materials; release archives do not.

## User Syntax

After installation, normal use is:

```stata
cicinference y x z, method(nosplit)
cicinference y x z, method(split)
cicinference y x z, method(kde)
cicinference y x z, method(bse) bootstrap(999) seed(2026)
cicinference y x z, method(bpc) bootstrap(999) seed(2026)
```

Results are stored in `r()` because the command is an inference/statistic
command rather than a full Stata estimator.
