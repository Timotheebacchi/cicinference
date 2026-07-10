version 16.0
clear all
set more off

capture confirm file "stata/ado/cicinference_bootstrap_common.ado"
if _rc {
    display as error "Run this test from the repository root."
    exit 601
}

capture confirm file "stata/tests/reference_outputs/reference_inputs.csv"
if _rc {
    display as error "Missing reference inputs. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

capture confirm file "stata/tests/reference_outputs/bootstrap_index_references.csv"
if _rc {
    display as error "Missing bootstrap index references. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

capture confirm file "stata/tests/reference_outputs/bootstrap_draw_references.csv"
if _rc {
    display as error "Missing bootstrap draw references. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

capture confirm file "stata/tests/reference_outputs/bootstrap_summary_references.csv"
if _rc {
    display as error "Missing bootstrap summary references. Run: Rscript stata/tests/make_reference_outputs.R"
    exit 601
}

adopath ++ "stata/ado"
do "stata/mata/cicinference_quantiles.mata"
do "stata/mata/cicinference_density.mata"
do "stata/mata/cicinference_variance.mata"
do "stata/mata/cicinference_bootstrap.mata"
do "stata/mata/cicinference_main.mata"

local inputs "stata/tests/reference_outputs/reference_inputs.csv"
local index_refs "stata/tests/reference_outputs/bootstrap_index_references.csv"
local draw_refs "stata/tests/reference_outputs/bootstrap_draw_references.csv"
local summary_refs "stata/tests/reference_outputs/bootstrap_summary_references.csv"
local tolerance 1e-9
local count_tolerance 1e-12

program define bootstrap_load_dataset
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

program define bootstrap_load_indices
    version 16.0
    syntax , Dataset(string) File(string)

    import delimited using "`file'", varnames(1) clear asdouble
    keep if dataset_id == "`dataset'"

    foreach role in Y X Z {
        preserve
        keep if role == "`role'"
        keep replicate draw index
        reshape wide index, i(replicate) j(draw)
        sort replicate
        ds index*
        local matname = lower("`role'") + "_idx"
        mkmat `r(varlist)', matrix(`matname')
        restore
    }
end

program define bootstrap_load_draw_reference
    version 16.0
    syntax , Dataset(string) File(string)

    import delimited using "`file'", varnames(1) clear asdouble
    keep if dataset_id == "`dataset'"
    sort replicate
    mkmat theta_boot, matrix(ref_boot_theta)
end

program define bootstrap_load_summary_scalars
    version 16.0
    syntax , Dataset(string) File(string)

    import delimited using "`file'", varnames(1) clear asdouble
    keep if dataset_id == "`dataset'"
    count
    if r(N) != 1 {
        display as error "Expected one bootstrap summary row for dataset=`dataset'."
        exit 459
    }

    scalar ref_theta_hat = theta_hat[1]
    scalar ref_se_boot = se_boot[1]
    scalar ref_bse_lower = bse_lower[1]
    scalar ref_bse_upper = bse_upper[1]
    scalar ref_bse_length = bse_length[1]
    scalar ref_bpc_lower = bpc_lower[1]
    scalar ref_bpc_upper = bpc_upper[1]
    scalar ref_bpc_length = bpc_length[1]
    scalar ref_B = B[1]
    scalar ref_n_y = n_y[1]
    scalar ref_n_x = n_x[1]
    scalar ref_n_z = n_z[1]
end

program define bootstrap_assert_scalar
    version 16.0
    syntax , Label(string) Observed(name) Expected(name) [ TOL(real 1e-9) ]

    scalar bootstrap_absdiff = abs(`observed' - `expected')
    if bootstrap_absdiff > `tol' {
        display as error "`label': observed = " %21.15g `observed' ///
            ", expected = " %21.15g `expected' ///
            ", abs diff = " %21.15g bootstrap_absdiff
        exit 9
    }
    display as text "`label': abs diff = " %9.2e bootstrap_absdiff
end

foreach dataset in xs_small xs_medium xs_edge_ties xs_edge_extreme {
    display as text "Checking common-index bootstrap for `dataset'"

    bootstrap_load_indices, dataset("`dataset'") file("`index_refs'")
    bootstrap_load_dataset, dataset("`dataset'") file("`inputs'")
    cicinference_bootstrap_common Y X Z, yindex(y_idx) xindex(x_idx) zindex(z_idx) level(95)

    matrix got_boot_theta = r(boot_theta)
    scalar got_theta_hat = r(theta_hat)
    scalar got_se_boot = r(se_boot)
    scalar got_bse_lower = r(bse_lower)
    scalar got_bse_upper = r(bse_upper)
    scalar got_bse_length = r(bse_length)
    scalar got_bpc_lower = r(bpc_lower)
    scalar got_bpc_upper = r(bpc_upper)
    scalar got_bpc_length = r(bpc_length)
    scalar got_B = r(B)
    scalar got_n_y = r(n_y)
    scalar got_n_x = r(n_x)
    scalar got_n_z = r(n_z)

    bootstrap_load_draw_reference, dataset("`dataset'") file("`draw_refs'")
    mata: st_numscalar("got_boot_max_absdiff", max(abs(st_matrix("got_boot_theta") :- st_matrix("ref_boot_theta"))))
    scalar zero_diff = 0
    bootstrap_assert_scalar, label("bootstrap draws `dataset'") observed(got_boot_max_absdiff) expected(zero_diff) tol(`tolerance')

    bootstrap_load_summary_scalars, dataset("`dataset'") file("`summary_refs'")
    bootstrap_assert_scalar, label("theta_hat `dataset'") observed(got_theta_hat) expected(ref_theta_hat) tol(`tolerance')
    bootstrap_assert_scalar, label("se_boot `dataset'") observed(got_se_boot) expected(ref_se_boot) tol(`tolerance')
    bootstrap_assert_scalar, label("bse_lower `dataset'") observed(got_bse_lower) expected(ref_bse_lower) tol(`tolerance')
    bootstrap_assert_scalar, label("bse_upper `dataset'") observed(got_bse_upper) expected(ref_bse_upper) tol(`tolerance')
    bootstrap_assert_scalar, label("bse_length `dataset'") observed(got_bse_length) expected(ref_bse_length) tol(`tolerance')
    bootstrap_assert_scalar, label("bpc_lower `dataset'") observed(got_bpc_lower) expected(ref_bpc_lower) tol(`tolerance')
    bootstrap_assert_scalar, label("bpc_upper `dataset'") observed(got_bpc_upper) expected(ref_bpc_upper) tol(`tolerance')
    bootstrap_assert_scalar, label("bpc_length `dataset'") observed(got_bpc_length) expected(ref_bpc_length) tol(`tolerance')
    bootstrap_assert_scalar, label("B `dataset'") observed(got_B) expected(ref_B) tol(`count_tolerance')
    bootstrap_assert_scalar, label("n_y `dataset'") observed(got_n_y) expected(ref_n_y) tol(`count_tolerance')
    bootstrap_assert_scalar, label("n_x `dataset'") observed(got_n_x) expected(ref_n_x) tol(`count_tolerance')
    bootstrap_assert_scalar, label("n_z `dataset'") observed(got_n_z) expected(ref_n_z) tol(`count_tolerance')
}

display as result "Common-index bootstrap reference comparisons completed."
