version 16.0
clear all
set more off

capture confirm file "stata/ado/cicinference.ado"
if _rc {
    display as error "Run this test from the repository root."
    exit 601
}

capture confirm file "stata/tests/reference_outputs/reference_inputs.csv"
if _rc {
    display as error "Missing reference inputs. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

capture confirm file "stata/tests/reference_outputs/reference_outputs.csv"
if _rc {
    display as error "Missing reference outputs. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

capture confirm file "stata/tests/reference_outputs/utility_nosplit_references.csv"
if _rc {
    display as error "Missing utility references. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

adopath ++ "stata/ado"

local inputs "stata/tests/reference_outputs/reference_inputs.csv"
local outputs "stata/tests/reference_outputs/reference_outputs.csv"
local utility_refs "stata/tests/reference_outputs/utility_nosplit_references.csv"
local tolerance 1e-9
local count_tolerance 1e-12

program define split_load_dataset
    version 16.0
    syntax , Dataset(string) File(string)

    import delimited using "`file'", varnames(1) clear asdouble
    keep if dataset_id == "`dataset'" & inlist(vector, "Y", "X", "Z")
    keep index vector value
    reshape wide value, i(index) j(vector) string
    rename valueY Y
    rename valueX X
    rename valueZ Z
    sort index
end

program define split_load_output_scalars
    version 16.0
    syntax , Dataset(string) File(string)

    import delimited using "`file'", varnames(1) clear asdouble
    keep if dataset_id == "`dataset'" & method == "split"
    count
    if r(N) != 1 {
        display as error "Expected one split output row for dataset=`dataset'."
        exit 459
    }

    scalar ref_theta_hat = theta_hat[1]
    scalar ref_se = se_implied[1]
    scalar ref_ci_lower = ci_lower[1]
    scalar ref_ci_upper = ci_upper[1]
    scalar ref_ci_length = ci_length[1]
    scalar ref_epsilon_n = epsilon_n[1]
    scalar ref_n_y = n_y[1]
    scalar ref_n_x = n_x[1]
    scalar ref_n_z = n_z[1]
end

program define split_load_utility_scalar
    version 16.0
    syntax , Dataset(string) Component(string) Scalar(name) File(string)

    import delimited using "`file'", varnames(1) clear asdouble
    keep if dataset_id == "`dataset'" & component == "`component'"
    count
    if r(N) != 1 {
        display as error "Expected one utility row for dataset=`dataset', component=`component'."
        exit 459
    }
    scalar `scalar' = value[1]
end

program define split_assert_scalar
    version 16.0
    syntax , Label(string) Observed(name) Expected(name) [ TOL(real 1e-9) ]

    scalar split_absdiff = abs(`observed' - `expected')
    if split_absdiff > `tol' {
        display as error "`label': observed = " %21.15g `observed' ///
            ", expected = " %21.15g `expected' ///
            ", abs diff = " %21.15g split_absdiff
        exit 9
    }
    display as text "`label': abs diff = " %9.2e split_absdiff
end

foreach dataset in xs_small xs_medium xs_edge_ties xs_edge_extreme {
    display as text "Checking split estimator for `dataset'"

    split_load_dataset, dataset("`dataset'") file("`inputs'")
    cicinference Y X Z, method(split) level(95)

    scalar got_theta_hat = r(theta_hat)
    scalar got_se = r(se)
    scalar got_ci_lower = r(ci_lower)
    scalar got_ci_upper = r(ci_upper)
    scalar got_ci_length = r(ci_length)
    scalar got_eps_hat = r(eps_hat)
    scalar got_eta = r(eta)
    scalar got_sigma_sq = r(sigma_sq)
    scalar got_epsilon_n = r(epsilon_n)
    scalar got_n_y = r(n_y)
    scalar got_n_x = r(n_x)
    scalar got_n_z = r(n_z)
    scalar got_n_half = r(n_half)

    split_load_output_scalars, dataset("`dataset'") file("`outputs'")
    split_load_utility_scalar, dataset("`dataset'") component("eps_hat") scalar(ref_eps_hat) file("`utility_refs'")
    split_load_utility_scalar, dataset("`dataset'") component("split_eta") scalar(ref_eta) file("`utility_refs'")
    split_load_utility_scalar, dataset("`dataset'") component("split_sigma_sq") scalar(ref_sigma_sq) file("`utility_refs'")
    split_load_utility_scalar, dataset("`dataset'") component("split_se") scalar(ref_se_utility) file("`utility_refs'")
    split_load_utility_scalar, dataset("`dataset'") component("split_n_half") scalar(ref_n_half) file("`utility_refs'")

    split_assert_scalar, label("theta_hat `dataset'") observed(got_theta_hat) expected(ref_theta_hat) tol(`tolerance')
    split_assert_scalar, label("se output `dataset'") observed(got_se) expected(ref_se) tol(`tolerance')
    split_assert_scalar, label("se utility `dataset'") observed(got_se) expected(ref_se_utility) tol(`tolerance')
    split_assert_scalar, label("ci_lower `dataset'") observed(got_ci_lower) expected(ref_ci_lower) tol(`tolerance')
    split_assert_scalar, label("ci_upper `dataset'") observed(got_ci_upper) expected(ref_ci_upper) tol(`tolerance')
    split_assert_scalar, label("ci_length `dataset'") observed(got_ci_length) expected(ref_ci_length) tol(`tolerance')
    split_assert_scalar, label("eps_hat `dataset'") observed(got_eps_hat) expected(ref_eps_hat) tol(`tolerance')
    split_assert_scalar, label("eta `dataset'") observed(got_eta) expected(ref_eta) tol(`tolerance')
    split_assert_scalar, label("sigma_sq `dataset'") observed(got_sigma_sq) expected(ref_sigma_sq) tol(`tolerance')
    split_assert_scalar, label("epsilon_n `dataset'") observed(got_epsilon_n) expected(ref_epsilon_n) tol(`tolerance')
    split_assert_scalar, label("n_y `dataset'") observed(got_n_y) expected(ref_n_y) tol(`count_tolerance')
    split_assert_scalar, label("n_x `dataset'") observed(got_n_x) expected(ref_n_x) tol(`count_tolerance')
    split_assert_scalar, label("n_z `dataset'") observed(got_n_z) expected(ref_n_z) tol(`count_tolerance')
    split_assert_scalar, label("n_half `dataset'") observed(got_n_half) expected(ref_n_half) tol(`count_tolerance')
}

display as result "Split estimator reference comparisons completed."
