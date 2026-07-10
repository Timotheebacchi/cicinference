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

capture confirm file "stata/tests/reference_outputs/reference_outputs.csv"
if _rc {
    display as error "Missing estimator references. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

do "stata/mata/cicinference_quantiles.mata"

local inputs "stata/tests/reference_outputs/reference_inputs.csv"
local utility_refs "stata/tests/reference_outputs/utility_nosplit_references.csv"
local outputs "stata/tests/reference_outputs/reference_outputs.csv"
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

program define phase3_load_nosplit_theta
    version 16.0
    syntax , Dataset(string) Mataname(name) File(string)

    import delimited using "`file'", varnames(1) clear asdouble
    keep if dataset_id == "`dataset'" & method == "no-split"
    count
    if r(N) != 1 {
        display as error "Expected one no-split reference row for dataset=`dataset'."
        exit 459
    }
    mata: `mataname' = st_data(., "theta_hat")
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
    display as text "Checking rank/quantile utilities for `dataset'"

    phase3_load_input_vector, dataset("`dataset'") vector("Y") mataname(y) file("`inputs'")
    phase3_load_input_vector, dataset("`dataset'") vector("X") mataname(x) file("`inputs'")
    phase3_load_input_vector, dataset("`dataset'") vector("Z") mataname(z) file("`inputs'")

    phase3_load_component, dataset("`dataset'") component("sorted_y") mataname(ref_sorted_y) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("sorted_x") mataname(ref_sorted_x) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("sorted_z") mataname(ref_sorted_z) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("uhat_ecdf_z_at_x") mataname(ref_uhat) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("left_quantile_index") mataname(ref_qindex) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("qcdf_transform") mataname(ref_qcdf) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("theta_hat_from_transform") mataname(ref_theta_transform) file("`utility_refs'")
    phase3_load_nosplit_theta, dataset("`dataset'") mataname(ref_theta_output) file("`outputs'")

    mata: computed_sorted_y = cic_sort_ascending(y)
    mata: computed_sorted_x = cic_sort_ascending(x)
    mata: computed_sorted_z = cic_sort_ascending(z)
    mata: computed_sorted_x_from_order = x[cic_order_ascending(x)]
    mata: computed_uhat = cic_ecdf_at(z, x)
    mata: computed_qindex = cic_left_quantile_indices(rows(y), ref_uhat)
    mata: computed_qcdf = cic_left_quantile(y, ref_uhat)
    mata: computed_transform = cic_rank_quantile_transform(y, x, z)
    mata: computed_transform_uhat = computed_transform[., 1]
    mata: computed_transform_qcdf = computed_transform[., 2]
    mata: computed_theta_transform = mean(computed_transform[., 2])
    mata: computed_theta_helper = cic_rank_quantile_mean(y, x, z)

    phase3_assert_close, label("sorted_y `dataset'") computed(computed_sorted_y) reference(ref_sorted_y) tol(`tolerance')
    phase3_assert_close, label("sorted_x `dataset'") computed(computed_sorted_x) reference(ref_sorted_x) tol(`tolerance')
    phase3_assert_close, label("sorted_z `dataset'") computed(computed_sorted_z) reference(ref_sorted_z) tol(`tolerance')
    phase3_assert_close, label("order_x `dataset'") computed(computed_sorted_x_from_order) reference(ref_sorted_x) tol(`tolerance')
    phase3_assert_close, label("ecdf_z_at_x `dataset'") computed(computed_uhat) reference(ref_uhat) tol(`tolerance')
    phase3_assert_close, label("left_quantile_index `dataset'") computed(computed_qindex) reference(ref_qindex) tol(`tolerance')
    phase3_assert_close, label("left_quantile_values `dataset'") computed(computed_qcdf) reference(ref_qcdf) tol(`tolerance')
    phase3_assert_close, label("transform_uhat `dataset'") computed(computed_transform_uhat) reference(ref_uhat) tol(`tolerance')
    phase3_assert_close, label("transform_qcdf `dataset'") computed(computed_transform_qcdf) reference(ref_qcdf) tol(`tolerance')
    phase3_assert_close, label("theta_transform `dataset'") computed(computed_theta_transform) reference(ref_theta_transform) tol(`tolerance')
    phase3_assert_close, label("theta_helper `dataset'") computed(computed_theta_helper) reference(ref_theta_output) tol(`tolerance')
}

display as result "Phase 3 rank/quantile utility tests passed."
