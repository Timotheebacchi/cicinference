{smcl}
{* *! version 0.1.0 10jul2026}{...}
{vieweralsosee "" "--"}{...}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col:{hi:cicinference} {hline 2}}Changes-in-Changes inference{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 17 2}
{cmd:cicinference} {it:y} {it:x} {it:z} {ifin}
[{cmd:,} {opt method(method)} {opt bootstrap(#)} {opt seed(#)} {opt level(#)}]

{pstd}
{it:method} is one of {cmd:nosplit}, {cmd:split}, {cmd:kde}, {cmd:bse}, or
{cmd:bpc}. The alias {cmd:no-split} is accepted for {cmd:nosplit}.

{title:Description}

{pstd}
{cmd:cicinference} estimates the Changes-in-Changes target
E[F_Y^{-1}(F_Z(X))] and reports a confidence interval using the selected
variance or bootstrap method from the R package {cmd:quantcdf.inference}.

{pstd}
The command uses complete observations of {it:y}, {it:x}, and {it:z};
Stata {cmd:if} and {cmd:in} qualifiers are respected. Current Stata support is
cross-sectional only.

{title:Options}

{phang}
{opt method(method)} selects the inference method. The default is
{cmd:method(nosplit)}.

{phang}
{opt bootstrap(#)} sets the number of bootstrap replications for {cmd:bse} and
{cmd:bpc}. Values below 200 are increased to 200 to match the R package.

{phang}
{opt seed(#)} sets Stata's RNG seed before the live Stata bootstrap. Matching an
R seed is not expected to produce identical bootstrap draws.

{phang}
{opt level(#)} sets the confidence level. Percent scale, such as {cmd:95}, and
fraction scale, such as {cmd:.95}, are both accepted.

{title:Stored results}

{pstd}
{cmd:cicinference} stores results in {cmd:r()} because it is an inference
statistic command rather than a full Stata estimation command.

{synoptset 20 tabbed}{...}
{synopt:{cmd:r(theta_hat)}}point estimate{p_end}
{synopt:{cmd:r(se)}}standard error, missing for percentile bootstrap{p_end}
{synopt:{cmd:r(ci_lower)}}lower confidence bound{p_end}
{synopt:{cmd:r(ci_upper)}}upper confidence bound{p_end}
{synopt:{cmd:r(ci_length)}}confidence interval length{p_end}
{synopt:{cmd:r(level)}}confidence level on fraction scale{p_end}
{synopt:{cmd:r(method)}}selected method{p_end}
{synopt:{cmd:r(n_y)}}Y sample size after missing-value handling{p_end}
{synopt:{cmd:r(n_x)}}X sample size after missing-value handling{p_end}
{synopt:{cmd:r(n_z)}}Z sample size after missing-value handling{p_end}
{synopt:{cmd:r(epsilon_n)}}bandwidth multiplier for nosplit, split, and kde{p_end}
{synopt:{cmd:r(se_boot)}}bootstrap standard deviation for bse and bpc{p_end}
{synopt:{cmd:r(B)}}bootstrap replications for bse and bpc{p_end}
{synopt:{cmd:r(result)}}method-specific result row vector{p_end}

{title:Examples}

{phang2}{cmd:. cicinference y x z, method(nosplit)}{p_end}
{phang2}{cmd:. cicinference y x z, method(split)}{p_end}
{phang2}{cmd:. cicinference y x z, method(kde)}{p_end}
{phang2}{cmd:. cicinference y x z, method(bse) bootstrap(999) seed(2026)}{p_end}
{phang2}{cmd:. cicinference y x z, method(bpc) bootstrap(999) seed(2026)}{p_end}

{title:Validation notes}

{pstd}
Development tests compare Mata/Stata outputs against deterministic R reference
CSVs in the source repository. Bootstrap validation uses common R-exported
resampling indices because R and Stata RNGs differ.

{pstd}
At the time this help file was written, Stata execution still had to be run on
a licensed local Stata installation before validation can be claimed.

{title:Reference}

{pstd}
Chhor, J., D'Haultfoeuille, X., L'Hour, J., and Mugnier, M. (2026).
Asymptotic Properties of Empirical Quantile-Based Estimators. Manuscript.
