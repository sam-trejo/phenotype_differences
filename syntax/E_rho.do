***CORRELATION TABLE

***********************************************************************************
*** LOAD WLS AND PGI DATA 
***********************************************************************************

matrix drop _all

***load data
use "${analytic_temp_plus_meta}", clear

***restrict variables
keep idpub is_grad pgi*

***hold data
tempfile hold
save `hold', replace

***********************************************************************************
*** CREATE LIST OF PHENOTYPES STRING
***********************************************************************************

global pgi_name ""
global pgi_num ""
global comma ""

local i=1

foreach var of varlist pgi* {
	
	label var `var' "`var'"
	display "`var'"
	global pgi_name "$pgi_name `var'"
 	global comma "$comma `var',"
	
	***add leading 0 to 1-9 for variable names
	local ilead = string(`i',"%02.0f")
	
	***rename polygenic scores
	rename `var' pgi`ilead'_
	
	global pgi_num "$pgi_num pgi`ilead'_"
	
	local i = `i' + 1
}

display "$pgi_name"
display "$pgi_num"
display "$comma"

global comma = substr("$comma", 1, length("${comma}")-2)

display "$comma"

***********************************************************************************
*** RESHAPE WIDE SO EACH ROW IS A SIBLING PAIR
***********************************************************************************

***reshape data so each row is a sibling pair
reshape wide pgi*, i(idpub) j(is_grad)

***********************************************************************************
*** RANDOMLY ASSIGN SIBLING ONE AND SIBLING TWO AND ESTIMATE CORRELATION 1000 TIMES
***********************************************************************************

***
gen y=.
gen x=.
foreach var in $pgi_num {
	
	***
	center `var'0 `var'1, standardize inplace
	
	***
	forval i=1/${corr_reps} { 
		quietly {
			***
			gen rand=rnormal()
			egen mn_rand=mean(rand), by(idpub)
			replace x=`var'0 if mn_rand<0
			replace y=`var'1 if mn_rand<0
			replace x=`var'1 if mn_rand>=0
			replace y=`var'0 if mn_rand>=0
			drop *rand
			
			***
			reg y x
			local b=_b[x]
			local se=(_se[x])^2
			matrix `var'=nullmat(`var') \ `b', `se' 
		}
	}
}

***********************************************************************************
*** COMBINES ESTIMATES
***********************************************************************************

***combine all phenotypes into one big matrix
foreach var in $pgi_num {
	mat combine = nullmat(combine), `var' 
}

***average sibling correlations within phenotype
mat row = J(rowsof(combine),1,1)
mat col = row'*combine
mat corr=col/rowsof(combine)

***convert variance to standard error
forvalues i=1/48 {
	mat out = nullmat(out), corr[1,2*`i'-1], corr[1,2*`i']^.5
}

***********************************************************************************
*** EXPORT CORRELATION DATA
***********************************************************************************

use `hold', clear

***convert rhos from matrix to variables
svmat corr

foreach var of varlist corr* {
	sum `var'
	replace `var' = `r(mean)' if `var'==.
}

local i=0
foreach var of varlist pgi* {
	
	***reattach phenotype names 
	local i = `i' + 1
	rename corr`i' rho_`var'
	
	local i = `i' + 1
	rename corr`i' var_rho_`var'
	gen se_rho_`var' = var_rho_`var'^.5
}

drop var*

save "${data}\clean\phen_diff_rho_${date}.dta", replace