version 16.0
clear all
set more off

capture confirm file "stata/mata/cicinference_quantiles.mata"
if _rc {
    display as error "Run this test from the repository root."
    exit 601
}

capture confirm file "stata/tests/reference_outputs/reference_inputs.csv"
if _rc {
    display as error "Missing reference inputs. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

capture confirm file "stata/tests/reference_outputs/utility_nosplit_references.csv"
if _rc {
    display as error "Missing utility references. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

do "stata/mata/cicinference_quantiles.mata"

local inputs "stata/tests/reference_outputs/reference_inputs.csv"
local utility_refs "stata/tests/reference_outputs/utility_nosplit_references.csv"
local tolerance 1e-11

program define phase3_load_input_vector
    version 16.0
    syntax , Dataset(string) Vector(string) Mataname(name) File(string)

    import delimited using "`file'", varnames(1) clear asdouble
    keep if dataset_id == "`dataset'" & vector == "`vector'"
    count
    if r(N) == 0 {
        display as error "No input rows for dataset=`dataset', vector=`vector'."
        exit 459
    }
    sort index
    mata: `mataname' = st_data(., "value")
end

program define phase3_load_component
    version 16.0
    syntax , Dataset(string) Component(string) Mataname(name) File(string)

    import delimited using "`file'", varnames(1) clear asdouble
    keep if dataset_id == "`dataset'" & component == "`component'"
    count
    if r(N) == 0 {
        display as error "No utility rows for dataset=`dataset', component=`component'."
        exit 459
    }
    sort index
    mata: `mataname' = st_data(., "value")
end

program define phase3_assert_close
    version 16.0
    syntax , Label(string) Computed(name) Reference(name) [ TOL(real 1e-11) ]

    mata: st_numscalar("phase3_same_shape", rows(`computed') == rows(`reference') & cols(`computed') == cols(`reference'))
    if phase3_same_shape != 1 {
        display as error "`label': dimension mismatch"
        exit 503
    }

    mata: st_numscalar("phase3_maxdiff", max(abs(`computed' :- `reference')))
    if phase3_maxdiff > `tol' {
        display as error "`label': max abs diff = " %21.15g phase3_maxdiff
        exit 9
    }
    display as text "`label': max abs diff = " %9.2e phase3_maxdiff
end

foreach dataset in xs_small xs_medium xs_edge_ties xs_edge_extreme {
    display as text "Checking deterministic no-split helpers for `dataset'"

    phase3_load_input_vector, dataset("`dataset'") vector("Y") mataname(y) file("`inputs'")
    phase3_load_component, dataset("`dataset'") component("fyhat_grid") mataname(ref_fyhat) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("ysortdiff") mataname(ref_ysortdiff) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("sorted_y") mataname(ref_sorted_y) file("`utility_refs'")

    mata: computed_fyhat = cic_empirical_grid(rows(y))
    mata: computed_ysortdiff = cic_adjacent_differences(y)
    mata: computed_ysortdiff_from_sorted = cic_adjacent_differences_sorted(ref_sorted_y)

    phase3_assert_close, label("fyhat_grid `dataset'") computed(computed_fyhat) reference(ref_fyhat) tol(`tolerance')
    phase3_assert_close, label("ysortdiff `dataset'") computed(computed_ysortdiff) reference(ref_ysortdiff) tol(`tolerance')
    phase3_assert_close, label("ysortdiff_sorted `dataset'") computed(computed_ysortdiff_from_sorted) reference(ref_ysortdiff) tol(`tolerance')
}

display as result "Phase 3 deterministic no-split helper tests passed."
