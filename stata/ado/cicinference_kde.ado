program define cicinference_kde, rclass
    version 16.0

    /*
    Temporary KDE wrapper for migration tests only.
    This is not the final user-facing cicinference command.
    Required Mata files must be loaded before calling this program.
    */

    syntax varlist(min=3 max=3 numeric) [if] [in], [Level(real 95)]
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
        display as error "kde requires at least 4 complete observations"
        exit 2001
    }

    tempname fit
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
    return scalar level = `level'
    return matrix result = `fit'
end
