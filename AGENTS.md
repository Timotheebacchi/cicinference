# AGENTS.md — Migration R package `cicinference` to Stata

## Role

You are a senior statistical software engineer specialized in R packages, Stata ado programming, Mata numerical programming, and econometric inference.

Your task is to migrate the existing R package `cicinference` into a Stata package with equivalent statistical behavior, equivalent user-facing functionality where possible, and carefully validated numerical outputs.

Do not blindly translate R line by line. Reconstruct the algorithms, identify their mathematical meaning, and implement them idiomatically in Stata/Mata.

## Main objective

Create a Stata implementation of the R package `cicinference`.

The final Stata package should expose a command such as:

```stata
cicinference y x z, method(nosplit)
cicinference y x z, method(split)
cicinference y x z, method(kde)
cicinference y x z, method(kde) bootstrap(999) seed(2026)
```

The exact syntax may be refined after inspecting the R package, but the Stata API must remain simple and natural for applied econometric users.

## Non-negotiable constraints

1. Preserve the statistical meaning of the R package.
2. Preserve numerical behavior as closely as possible.
3. Do not implement expensive numerical routines in slow ado loops if they should be in Mata.
4. Use Mata for sorting, quantile computation, two-pointer search, cumulative sums, bootstrap internals, variance computation, and density or KDE-like operations.
5. Use ado files only for Stata user interface, option parsing, data extraction, returned results, and display.
6. Do not introduce new statistical assumptions unless explicitly documented.
7. Do not silently change defaults from the R package.
8. Do not remove functionality merely because it is inconvenient to port.
9. If exact equivalence is impossible, document the discrepancy precisely.
10. Every migrated component must have a test comparing it against the R reference output.

## Expected repository structure

Create a Stata-side structure similar to:

```text
stata/
  ado/
    cicinference.ado
    cicinference.sthlp

  mata/
    cicinference_utils.mata
    cicinference_quantiles.mata
    cicinference_density.mata
    cicinference_variance.mata
    cicinference_bootstrap.mata
    cicinference_main.mata

  tests/
    make_reference_outputs.R
    test_nosplit.do
    test_split.do
    test_kde.do
    test_bootstrap.do
    benchmark.do

  examples/
    basic_usage.do
    replication_example.do

  README-stata.md
```

If Stata packaging conventions require a different layout, propose the change before applying it.

## Work plan

### Phase 1 — Inspect and map the R package

Before writing Stata code, inspect the entire repository.

Produce a file:

```text
stata/MIGRATION_MAP.md
```

It must contain:

1. List of exported R functions.
2. List of internal R functions.
3. Function dependency graph.
4. Identification of pure algorithmic functions.
5. Identification of interface/display/documentation functions.
6. Identification of C++/Rcpp functions, if any.
7. For each major function:

   * inputs;
   * outputs;
   * side effects;
   * mathematical role;
   * computational complexity;
   * proposed Stata/Mata target file.

Do not start implementation before this map exists.

### Phase 2 — Create R reference outputs

Create deterministic reference tests from the R package.

Add:

```text
stata/tests/make_reference_outputs.R
```

This script must:

1. Load the R package locally.
2. Generate small, medium, and moderately large datasets.
3. Run all main methods.
4. Save reference outputs in a simple exchange format readable by Stata, preferably CSV or JSON.
5. Include fixed seeds.
6. Include edge cases:

   * small n;
   * ties in Y, X, or Z;
   * missing values if the R package handles them;
   * extreme quantile values;
   * bootstrap with small B;
   * KDE method if available.

The Stata implementation will be judged against these reference outputs.

### Phase 3 — Implement low-level Mata utilities

Implement and test utilities first.

Priority order:

1. Sorting helpers.
2. Left-continuous quantile function.
3. Empirical CDF/rank helpers.
4. Two-pointer interval count routine.
5. Cumulative sum helpers.
6. Kernel/density helpers.
7. Bootstrap sampling helpers.
8. Any influence-function or variance helper.

Each utility must have a small test.

Avoid cleverness unless it improves asymptotic complexity or numerical stability.

### Phase 4 — Implement main estimators

Implement methods in this order:

1. `nosplit`
2. `split`
3. `kde`
4. bootstrap variants

For each method:

1. First implement the simplest correct version.
2. Compare against R reference outputs.
3. Only then optimize.
4. Keep the optimized version readable.
5. Add comments explaining the mathematical object being computed.

### Phase 5 — Implement ado interface

Create:

```text
stata/ado/cicinference.ado
```

The ado command must:

1. Parse variables and options.
2. Respect Stata `if` and `in` qualifiers.
3. Respect missing-value handling.
4. Pass clean numeric vectors to Mata.
5. Display results in a standard Stata format.
6. Store results in `r()` or `e()` as appropriate.

Prefer `e()` if the command behaves like an estimator with standard errors and confidence intervals. Prefer `r()` if it behaves like a test/statistic command. Decide explicitly and document the choice.

Potential options:

```stata
method(nosplit|split|kde)
bootstrap(integer)
seed(integer)
level(real)
timings
```

Add more options only if they exist in the R package or are necessary for equivalence.

### Phase 6 — Documentation

Create:

```text
stata/ado/cicinference.sthlp
stata/README-stata.md
```

The help file must include:

1. Syntax.
2. Description.
3. Options.
4. Stored results.
5. Examples.
6. Notes on equivalence with the R package.
7. References if present in the R package.

The README must explain how to install and test the Stata implementation locally.

### Phase 7 — Validation

Create Stata test scripts that compare outputs against R references.

Required tests:

```text
stata/tests/test_nosplit.do
stata/tests/test_split.do
stata/tests/test_kde.do
stata/tests/test_bootstrap.do
stata/tests/benchmark.do
```

Each test must:

1. Load reference data.
2. Run the Stata command.
3. Compare estimates, standard errors, confidence intervals, p-values, or relevant internal quantities.
4. Use tolerances appropriate for floating-point arithmetic.
5. Fail loudly if discrepancies exceed tolerance.

Use strict tolerances for deterministic non-bootstrap quantities. Use looser tolerances only when justified.

### Phase 8 — Performance

After correctness is established, benchmark the Stata version.

Measure:

1. Runtime for n = 1,000.
2. Runtime for n = 10,000.
3. Runtime for n = 100,000 if feasible.
4. Runtime by method.
5. Bootstrap runtime for small and moderate B.

Document complexity in:

```text
stata/PERFORMANCE_NOTES.md
```

Do not sacrifice correctness for performance.

## Mathematical caution

Be especially careful with:

1. Indexing: R is 1-indexed; Mata is also 1-indexed, but Stata data handling differs.
2. Quantile conventions: preserve the package’s left-continuous quantile convention if used.
3. Ties: do not assume continuous data unless the R package assumes it.
4. Missing values: Stata missing values are larger than all real numbers; handle them explicitly.
5. Random seeds: R and Stata RNGs differ, so bootstrap equality may require comparing distributions or using exported bootstrap indices.
6. Floating-point tolerance: exact equality is usually inappropriate for doubles.
7. Sorting stability: check whether the R implementation relies on stable sorting.
8. KDE/bandwidth rules: preserve exactly if possible.

## Implementation standards

Use clear names.

Bad:

```mata
a = b[.,1]
```

Better:

```mata
sorted_y = data[order_index, 1]
```

Use comments to explain nontrivial formulas.

Do not over-comment trivial assignments.

Prefer small functions with explicit contracts.

Every Mata function should say:

```text
Inputs:
Outputs:
Assumptions:
```

## Git discipline

Before making large changes:

1. Check current git status.
2. Do not overwrite user work.
3. Prefer small commits or at least small coherent diffs.
4. After each phase, summarize changed files and remaining issues.

Do not reformat unrelated R files.

Do not rename existing R functions unless necessary.

## Definition of done

The migration is done only when:

1. The Stata command runs on example datasets.
2. All main R methods have Stata equivalents.
3. Reference tests pass within documented tolerances.
4. Help file exists and is usable from Stata.
5. Stored results are documented.
6. Performance notes exist.
7. Known differences from R are explicitly listed.
8. The final diff is reviewed for numerical, statistical, and Stata-interface errors.

## First task

Start by inspecting the repository.

Do not write implementation code yet.

Create `stata/MIGRATION_MAP.md` containing the architecture map and migration plan described above. Then stop and report:

1. What the R package does.
2. Which functions are hardest to migrate.
3. Which functions should be implemented in Mata.
4. Which functions should remain ado-level.
5. A proposed first implementation milestone.
