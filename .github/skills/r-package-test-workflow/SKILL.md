---
name: r-package-test-workflow
description: "Use when you need to verify R package changes with testthat. Create or update focused regression tests, run the relevant test file or package test suite, and interpret the result to confirm the change works."
argument-hint: "Target file, function, or behavior to verify"
user-invocable: true
---

# R Package Test Workflow

Use this skill when you need to confirm that a change in an R package behaves correctly by adding or running `testthat` tests.

## When to Use
- A function was changed and needs a regression test.
- You want to confirm a bug fix or feature behaves as expected.
- You need a small, targeted test instead of a full package-wide rewrite.

## Procedure
1. Identify the changed behavior, the public function or helper it affects, and the smallest reproducible input.
2. Check the existing test structure in `tests/testthat/` and reuse the current style, helpers, and expectations.
3. Add or update a focused `testthat::test_that()` case that proves the behavior before and after the change.
4. Prefer explicit assertions about output structure, values, warnings, and errors over broad snapshot checks.
5. Run the narrowest useful validation first: a single test file, then the package test suite if needed.
6. If the test fails, decide whether the code or the expectation is wrong, then update only the local slice and rerun the same check.
7. Confirm the final result with a passing test run and keep the test readable enough to serve as future documentation.

## Quality Checks
- The test fails for the original bug or missing behavior.
- The test passes after the fix.
- The assertion is specific enough to catch regressions without being fragile.
- The test matches the repo's current `testthat` conventions.

## References
- Existing package tests: `./tests/testthat/test-cic.R`
- Package metadata: `./DESCRIPTION`
- CI test coverage workflow: `./.github/workflows/test-coverage.yaml`