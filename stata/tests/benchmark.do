version 16.0
clear all
set more off

capture confirm file "stata/ado/cicinference.ado"
if _rc {
    display as error "Run this benchmark from the repository root."
    exit 601
}

adopath ++ "stata/ado"

program define benchmark_one_size
    version 16.0
    syntax , N(integer)

    clear
    set obs `n'
    set seed 2026
    generate double Y = rnormal()
    generate double X = rnormal()
    generate double Z = rnormal()

    foreach method in nosplit split kde {
        timer clear
        timer on 1
        quietly cicinference Y X Z, method(`method') level(95)
        timer off 1
        display as text "n=`n', method=`method'"
        timer list 1
    }

    timer clear
    timer on 1
    quietly cicinference Y X Z, method(bse) bootstrap(200) seed(2026) level(95)
    timer off 1
    display as text "n=`n', method=bse, B=200"
    timer list 1

    timer clear
    timer on 1
    quietly cicinference Y X Z, method(bpc) bootstrap(200) seed(2026) level(95)
    timer off 1
    display as text "n=`n', method=bpc, B=200"
    timer list 1
end

benchmark_one_size, n(1000)
benchmark_one_size, n(10000)

display as result "cicinference Stata benchmark completed."
