 
***********************************************************************************
*** PREP DATA
***********************************************************************************

use "${analytic}", clear

 robvar meta_pgi_alive75 if dna, by(dead)
 robvar meta_pgi_span if dna, by(dead) 
 
 local i=0
foreach var of global stub2 {
	quietly robvar pgi_`var' if dna, by(dead) 
	local p=round(`r(p_w50)', .01)
	display "`var': p=`p'"
 }
 

***********************************************************************************
*** PREP DATA
***********************************************************************************

***
global vars out_bmi out_hgt out_cog out_edu out_extra out_neuro out_open
keep if dna_sib==0

***********************************************************************************
*** SD RATIO PROGRAM
***********************************************************************************

program define sd_ratio, rclass
	args v1 v2
	confirm variable `v1'
	confirm variable `v2'
	
	quietly sum `v1' if `v2' == 0 
	local num = `r(Var)'
	
	quietly sum `v1' if `v2' == 1 
	local den =  `r(Var)'

	return scalar ratio = (`num'/`den')
	display (`num'/`den')
end

***********************************************************************************
*** COLUMNS 1 & 2
***********************************************************************************

preserve
	keep if dna==0
	foreach var of varlist $vars {
		sum `var'

		replace `var' = (`r(Var)')

	}
	eststo: estpost sum	$vars
restore

***
preserve
	keep if dna==0
	foreach var of varlist $vars {
		sum `var'
		replace `var' = `r(N)'
	}
	eststo: estpost sum	$vars
restore

***********************************************************************************
*** COLUMNS 3 & 4
***********************************************************************************

preserve
	keep if dna==1
	foreach var of varlist $vars {
		sum `var'
		replace `var' = (`r(Var)')
	}
	eststo: estpost sum	$vars
restore

***
preserve
	keep if dna==1
	foreach var of varlist $vars {
		sum `var'
		replace `var' = `r(N)'
	}	
	eststo: estpost sum	$vars
restore

***********************************************************************************
*** COLUMN 5
***********************************************************************************

preserve
	foreach var of varlist $vars {
		sd_ratio `var' dna
		replace `var' = `r(ratio)'
	}
	eststo: estpost sum	$vars
restore

***********************************************************************************
*** COLUMN 6
***********************************************************************************

preserve
	foreach var of varlist $vars {
		display "`var'"
		robvar `var', by(dna)
		replace `var' = `r(p_w50)'
	}
	eststo: estpost sum	$vars
restore

***********************************************************************************
*** 
***********************************************************************************	
		
esttab using "{$tables}_", /// "${table}\var_test_${date}.tex",
			 cells("mean(fmt(a2))") ///
			 noobs nonumber nodepvar label replace booktabs fragment ///
			 collabels(none) ///
			 mlabels("$Var_1$" "$N_1$" "$Var_2$" "$N_2$" "$\frac{Var_2}{Var_1}$" " ") ///
			 mgroups("\textbf{Not Genotyped}" "\textbf{Genotyped}" "\textbf{Ratio}" "\textbf{\textit{p}-value}", pattern(1 0 1 0 1 1) span prefix(\multicolumn{@span}{c}{) suffix(}))
