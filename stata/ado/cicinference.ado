program define cicinference, rclass
    version 16.0

    syntax varlist(min=3 max=3 numeric) [if] [in], ///
        [Method(string) Bootstrap(integer 1000) Seed(integer) Level(real 95)]
    marksample touse
    markout `touse' `varlist'

    if "`method'" == "" {
        local method "nosplit"
    }
    local method_norm = lower(strtrim("`method'"))
    local method_norm : subinstr local method_norm "-" "", all

    if !inlist("`method_norm'", "nosplit", "split", "kde", "bse", "bpc") {
        display as error "method() must be one of nosplit, split, kde, bse, or bpc"
        exit 198
    }

    if (`level' > 1) {
        local level = `level' / 100
    }
    if (`level' <= 0 | `level' >= 1) {
        display as error "level() must be in (0,1) or a percentage in (0,100)"
        exit 198
    }

    quietly count if `touse'
    if (r(N) < 4) {
        display as error "cicinference requires at least 4 complete observations"
        exit 2001
    }

    tempname loadcheck
    capture mata: st_numscalar("`loadcheck'", ///
        cic_min3(1, 2, 3) + ///
        cic_bootstrap_summary_from_values(0, (1 \ 2), 0.95)[1,1])
    if _rc {
        findfile "cicinference.ado"
        local cic_ado_path "`r(fn)'"
        local cic_ado_dir = substr("`cic_ado_path'", 1, strrpos("`cic_ado_path'", "/") - 1)
        local cic_stata_dir = substr("`cic_ado_dir'", 1, strrpos("`cic_ado_dir'", "/") - 1)

        foreach matafile in ///
            cicinference_quantiles.mata ///
            cicinference_density.mata ///
            cicinference_variance.mata ///
            cicinference_bootstrap.mata ///
            cicinference_main.mata {
            capture findfile "`matafile'"
            if _rc {
                local mata_path "`cic_stata_dir'/mata/`matafile'"
                capture confirm file "`mata_path'"
                if _rc {
                    display as error "could not find required Mata file `matafile'"
                    exit 601
                }
            }
            else {
                local mata_path "`r(fn)'"
            }
            quietly do "`mata_path'"
        }
    }

    gettoken y rest : varlist
    gettoken x z : rest

    tempname fit
    local method_label "`method_norm'"
    local display_se "."
    local epsilon_n .
    local B_used .

    if "`method_norm'" == "nosplit" {
        local method_label "no-split"
        mata: st_matrix("`fit'", cic_nosplit_fit( ///
            st_data(., "`y'", "`touse'"), ///
            st_data(., "`x'", "`touse'"), ///
            st_data(., "`z'", "`touse'"), ///
            `level'))
        matrix colnames `fit' = theta_hat se ci_lower ci_upper ci_length ///
            eps_hat eta sigma_sq epsilon_n n_y n_x n_z

        return scalar theta_hat = `fit'[1,1]
        return scalar se = `fit'[1,2]
        return scalar ci_lower = `fit'[1,3]
        return scalar ci_upper = `fit'[1,4]
        return scalar ci_length = `fit'[1,5]
        return scalar eps_hat = `fit'[1,6]
        return scalar eta = `fit'[1,7]
        return scalar sigma_sq = `fit'[1,8]
        return scalar epsilon_n = `fit'[1,9]
        return scalar n_y = `fit'[1,10]
        return scalar n_x = `fit'[1,11]
        return scalar n_z = `fit'[1,12]
    }
    else if "`method_norm'" == "split" {
        local method_label "split"
        mata: st_matrix("`fit'", cic_split_fit( ///
            st_data(., "`y'", "`touse'"), ///
            st_data(., "`x'", "`touse'"), ///
            st_data(., "`z'", "`touse'"), ///
            `level'))
        matrix colnames `fit' = theta_hat se ci_lower ci_upper ci_length ///
            eps_hat eta sigma_sq epsilon_n n_y n_x n_z n_half

        return scalar theta_hat = `fit'[1,1]
        return scalar se = `fit'[1,2]
        return scalar ci_lower = `fit'[1,3]
        return scalar ci_upper = `fit'[1,4]
        return scalar ci_length = `fit'[1,5]
        return scalar eps_hat = `fit'[1,6]
        return scalar eta = `fit'[1,7]
        return scalar sigma_sq = `fit'[1,8]
        return scalar epsilon_n = `fit'[1,9]
        return scalar n_y = `fit'[1,10]
        return scalar n_x = `fit'[1,11]
        return scalar n_z = `fit'[1,12]
        return scalar n_half = `fit'[1,13]
    }
    else if "`method_norm'" == "kde" {
        local method_label "kde"
        mata: st_matrix("`fit'", cic_kde_fit( ///
            st_data(., "`y'", "`touse'"), ///
            st_data(., "`x'", "`touse'"), ///
            st_data(., "`z'", "`touse'"), ///
            `level'))
        matrix colnames `fit' = theta_hat se ci_lower ci_upper ci_length ///
            eps_hat eta sigma_sq epsilon_n n_y n_x n_z

        return scalar theta_hat = `fit'[1,1]
        return scalar se = `fit'[1,2]
        return scalar ci_lower = `fit'[1,3]
        return scalar ci_upper = `fit'[1,4]
        return scalar ci_length = `fit'[1,5]
        return scalar eps_hat = `fit'[1,6]
        return scalar eta = `fit'[1,7]
        return scalar sigma_sq = `fit'[1,8]
        return scalar epsilon_n = `fit'[1,9]
        return scalar n_y = `fit'[1,10]
        return scalar n_x = `fit'[1,11]
        return scalar n_z = `fit'[1,12]
    }
    else {
        local B_used = `bootstrap'
        if (`B_used' < 200) {
            display as text "warning: B < 200; using B = 200 to match the R package"
            local B_used 200
        }
        if "`seed'" != "" {
            set seed `seed'
            return scalar seed = `seed'
        }

        mata: st_matrix("`fit'", cic_bootstrap_fit( ///
            st_data(., "`y'", "`touse'"), ///
            st_data(., "`x'", "`touse'"), ///
            st_data(., "`z'", "`touse'"), ///
            `B_used', ///
            `level'))
        matrix colnames `fit' = theta_hat se_boot bse_lower bse_upper bse_length ///
            bpc_lower bpc_upper bpc_length B n_y n_x n_z

        return scalar theta_hat = `fit'[1,1]
        return scalar se_boot = `fit'[1,2]
        return scalar B = `fit'[1,9]
        return scalar n_y = `fit'[1,10]
        return scalar n_x = `fit'[1,11]
        return scalar n_z = `fit'[1,12]

        if "`method_norm'" == "bse" {
            local method_label "bse"
            return scalar se = `fit'[1,2]
            return scalar ci_lower = `fit'[1,3]
            return scalar ci_upper = `fit'[1,4]
            return scalar ci_length = `fit'[1,5]
        }
        else {
            local method_label "bpc"
            return scalar se = .
            return scalar ci_lower = `fit'[1,6]
            return scalar ci_upper = `fit'[1,7]
            return scalar ci_length = `fit'[1,8]
        }
    }

    return scalar level = `level'
    return local method "`method_label'"
    return matrix result = `fit'

    display as text _newline "Changes-in-Changes inference"
    display as text "Method: " as result "`method_label'"
    display as text "Complete observations: " as result %9.0g r(n_x)
    if inlist("`method_norm'", "bse", "bpc") {
        display as text "Bootstrap replications: " as result %9.0g r(B)
    }

    display as text _newline _col(3) "theta_hat" ///
        _col(18) "std. err." ///
        _col(33) "ci lower" ///
        _col(48) "ci upper" ///
        _col(63) "ci length"
    display as result _col(3) %10.6g r(theta_hat) ///
        _col(18) %10.6g r(se) ///
        _col(33) %10.6g r(ci_lower) ///
        _col(48) %10.6g r(ci_upper) ///
        _col(63) %10.6g r(ci_length)
end
