version 16.0
clear all
set more off

capture confirm file "stata/mata/cicinference_quantiles.mata"
if _rc {
    display as error "Run this test from the repository root."
    exit 601
}

capture confirm file "stata/mata/cicinference_density.mata"
if _rc {
    display as error "Missing density Mata utilities."
    exit 601
}

capture confirm file "stata/mata/cicinference_variance.mata"
if _rc {
    display as error "Missing variance Mata utilities."
    exit 601
}

capture confirm file "stata/tests/reference_outputs/utility_nosplit_references.csv"
if _rc {
    display as error "Missing utility references. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

do "stata/mata/cicinference_quantiles.mata"
do "stata/mata/cicinference_density.mata"
do "stata/mata/cicinference_variance.mata"

local utility_refs "stata/tests/reference_outputs/utility_nosplit_references.csv"
local tolerance 1e-10

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
    syntax , Label(string) Computed(name) Reference(name) [ TOL(real 1e-10) ]

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
    display as text "Checking Phase 3-b density/variance helpers for `dataset'"

    phase3_load_component, dataset("`dataset'") component("sorted_uhat") mataname(ref_sorted_uhat) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("fyhat_grid") mataname(ref_fyhat) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("epsilon_n_default") mataname(ref_epsilon) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("density_bandwidth") mataname(ref_bandwidth) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("rect_counts") mataname(ref_counts) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("fuhat_density") mataname(ref_density) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("ysortdiff") mataname(ref_ysortdiff) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("fast_eta_u") mataname(ref_fast_eta_u) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("fast_eta_reverse_cumsum") mataname(ref_fast_eta_reverse_cumsum) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("fast_eta_delta_v") mataname(ref_fast_eta_delta_v) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("fast_eta_term2") mataname(ref_fast_eta_term2) file("`utility_refs'")
    phase3_load_component, dataset("`dataset'") component("fast_eta_self") mataname(ref_fast_eta_self) file("`utility_refs'")

    mata: epsilon_n = ref_epsilon[1]
    mata: computed_bandwidth = cic_density_bandwidth(ref_fyhat, epsilon_n)
    mata: computed_counts = cic_rect_counts_sorted(ref_sorted_uhat, ref_fyhat, ref_bandwidth)
    mata: computed_density_from_counts = cic_counts_to_density(ref_counts, ref_bandwidth, rows(ref_sorted_uhat))
    mata: computed_density_from_bandwidths = cic_rect_density_from_bandwidths(ref_sorted_uhat, ref_fyhat, ref_bandwidth)
    mata: computed_density_from_epsilon = cic_rect_density_sorted(ref_sorted_uhat, ref_fyhat, epsilon_n)
    mata: computed_fast_eta_u = ref_ysortdiff :* ref_density
    mata: computed_fast_eta_reverse_cumsum = cic_reverse_cumsum(ref_fast_eta_u)
    mata: computed_fast_eta_delta_v = cic_grid_increments(ref_fyhat)
    mata: computed_fast_eta_term2 = sum(ref_fast_eta_u :* ref_fyhat) * sum(ref_fast_eta_u :* ref_fyhat)
    mata: computed_fast_eta = cic_fast_eta(ref_ysortdiff, ref_density, ref_ysortdiff, ref_density, ref_fyhat)
    mata: computed_fast_eta_self = cic_fast_eta_self(ref_ysortdiff, ref_density, ref_fyhat)

    phase3_assert_close, label("density_bandwidth `dataset'") computed(computed_bandwidth) reference(ref_bandwidth) tol(`tolerance')
    phase3_assert_close, label("rect_counts `dataset'") computed(computed_counts) reference(ref_counts) tol(`tolerance')
    phase3_assert_close, label("density_from_counts `dataset'") computed(computed_density_from_counts) reference(ref_density) tol(`tolerance')
    phase3_assert_close, label("density_from_bandwidths `dataset'") computed(computed_density_from_bandwidths) reference(ref_density) tol(`tolerance')
    phase3_assert_close, label("density_from_epsilon `dataset'") computed(computed_density_from_epsilon) reference(ref_density) tol(`tolerance')
    phase3_assert_close, label("fast_eta_u `dataset'") computed(computed_fast_eta_u) reference(ref_fast_eta_u) tol(`tolerance')
    phase3_assert_close, label("fast_eta_reverse_cumsum `dataset'") computed(computed_fast_eta_reverse_cumsum) reference(ref_fast_eta_reverse_cumsum) tol(`tolerance')
    phase3_assert_close, label("fast_eta_delta_v `dataset'") computed(computed_fast_eta_delta_v) reference(ref_fast_eta_delta_v) tol(`tolerance')
    phase3_assert_close, label("fast_eta_term2 `dataset'") computed(computed_fast_eta_term2) reference(ref_fast_eta_term2) tol(`tolerance')
    phase3_assert_close, label("fast_eta `dataset'") computed(computed_fast_eta) reference(ref_fast_eta_self) tol(`tolerance')
    phase3_assert_close, label("fast_eta_self `dataset'") computed(computed_fast_eta_self) reference(ref_fast_eta_self) tol(`tolerance')
}

display as result "Phase 3-b density/variance helper tests completed."
