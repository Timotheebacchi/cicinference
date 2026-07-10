program define cicinference_bootstrap_common, rclass
    version 16.0

    /*
    Temporary bootstrap wrapper for migration tests only.
    It uses common R-exported bootstrap index matrices, so results are
    deterministic and comparable despite R/Stata RNG differences.
    */

    syntax varlist(min=3 max=3 numeric) [if] [in], ///
        YIndex(name) XIndex(name) ZIndex(name) [Level(real 95)]
    marksample touse
    markout `touse' `varlist'

    gettoken y rest : varlist
    gettoken x z : rest

    if (`level' > 1) {
        local level = `level' / 100
    }
    if (`level' <= 0 | `level' >= 1) {
        display as error "level() must be in (0,1) or a percentage in (0,100)"
        exit 198
    }

    quietly count if `touse'
    if (r(N) < 4) {
        display as error "bootstrap requires at least 4 complete observations"
        exit 2001
    }

    tempname fit boot
    mata: st_matrix("`fit'", cic_bootstrap_fit_with_indices( ///
        st_data(., "`y'", "`touse'"), ///
        st_data(., "`x'", "`touse'"), ///
        st_data(., "`z'", "`touse'"), ///
        st_matrix("`yindex'"), ///
        st_matrix("`xindex'"), ///
        st_matrix("`zindex'"), ///
        `level'))
    mata: st_matrix("`boot'", cic_bootstrap_values_with_indices( ///
        cic_sort_ascending(st_data(., "`y'", "`touse'")), ///
        st_data(., "`x'", "`touse'"), ///
        cic_sort_ascending(st_data(., "`z'", "`touse'")), ///
        st_matrix("`yindex'"), ///
        st_matrix("`xindex'"), ///
        st_matrix("`zindex'")))

    matrix colnames `fit' = theta_hat se_boot bse_lower bse_upper bse_length ///
        bpc_lower bpc_upper bpc_length B n_y n_x n_z

    return scalar theta_hat = `fit'[1,1]
    return scalar se_boot = `fit'[1,2]
    return scalar bse_lower = `fit'[1,3]
    return scalar bse_upper = `fit'[1,4]
    return scalar bse_length = `fit'[1,5]
    return scalar bpc_lower = `fit'[1,6]
    return scalar bpc_upper = `fit'[1,7]
    return scalar bpc_length = `fit'[1,8]
    return scalar B = `fit'[1,9]
    return scalar n_y = `fit'[1,10]
    return scalar n_x = `fit'[1,11]
    return scalar n_z = `fit'[1,12]
    return scalar level = `level'
    return matrix result = `fit'
    return matrix boot_theta = `boot'
end
