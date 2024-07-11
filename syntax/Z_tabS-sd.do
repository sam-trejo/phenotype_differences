 
***********************************************************************************
*** SD TESTS ON GENOTYPED SAMPLE WHO ARE DEAD/ALIVE by 2020
***********************************************************************************

use "${analytic}", clear

 robvar pgi_meta if dna, by(dead)
 
 local i=0
foreach var of global stub2 {
	quietly robvar pgi_`var' if dna, by(dead) 
	local p=round(`r(p_w50)', .01)
	display "`var': p=`p'"
 }
 

***********************************************************************************
*** PREP DATA TO EXPORT PHENOTYPIC SD TEST FOR TABLE IN SI
***********************************************************************************

estimates clear 

***
keep if dna_sib==0

***
global vars out_bmi out_hgt out_cog out_edu out_extra out_neuro out_open

foreach var of global vars {
	egen count = count(`var'), by(idpub)	
	replace `var' = . if count!=2
	drop count
}

***********************************************************************************
*** SD RATIO PROGRAM
***********************************************************************************


program drop _all
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
*** COLUMNS 1 
***********************************************************************************

preserve
	keep if dna==0
	foreach var of varlist $vars {
		sum `var'

		replace `var' = (`r(Var)')

	}
	eststo: estpost sum	$vars
restore

***********************************************************************************
*** COLUMNS 2 & 3
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
*** COLUMN 4
***********************************************************************************

preserve
	foreach var of varlist $vars {
		sd_ratio `var' dna
		replace `var' = `r(ratio)'
	}
	eststo: estpost sum	$vars
restore

***********************************************************************************
*** COLUMN 5
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
		
esttab using "${table}\tabS-sd_${date}.tex", ///
			 cells("mean(fmt(a2))") ///
			 noobs nonumber nodepvar label replace booktabs fragment ///
			 collabels(none) ///
			 mlabels("$var(y_{2j})$" "$var(y_{1j})$" "$N$ Sibling Pairs" "$\frac{var(y_{2j})}{var(y_{1j})}$" "\textit{p}-value")
			 