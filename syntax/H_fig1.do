***********************************************************************************
*** LOAD PROGRAMS
***********************************************************************************

do "${syntax}/programs/phendiff.do"
do "${syntax}/programs/bky2006.do"

***********************************************************************************
*** PANEL B: ONE RANDOM SIBLING SAMPLE (SINGLE REP)
***********************************************************************************

use "${analytic}", clear

matrix drop _all

preserve 
	***randomly select a focal sib
	gen rand = rnormal() if dna_sib==1
	egen fam_mn_rand = mean(rand), by(idpub)
	gen focal = rand>fam_mn_rand
	drop rand fam_mn_rand
		
	global rowname ""

	foreach var of global stub {
		display "`var'"
		
		global rowname "$rowname `var'"

		quietly reghdfe out_`var' pgi_`var' ///
						if dna_sib==1 ///
						, absorb(idpub)
		local b_fe = _b[pgi_`var']
		local se_fe = _se[pgi_`var']
		
		replace pgi_`var'=. if focal==0		
		quietly sum rho_pgi_`var'
		quietly phendiff out_`var' pgi_`var' ///
				 if dna_sib==1 ///
				 , fam(idpub) ///
				 rho(`r(mean)') ///
				 delta(2104)
		
		mat betas = nullmat(betas) \ `b_fe', `se_fe', `r(beta)', `r(se)'
	}			

	mat rownames betas = $rowname			
	mat colnames betas = "b_fe" "se_fe" "b_pd" "se_pd"
		
	***
	svmat2 betas, names(col) rnames(pheno) 

	keep b_fe se_fe b_pd se_pd pheno
	keep if pheno!=""
	
	***generate confidence interval variables
	gen b_fe_ci_hi = b_fe + 1.96*se_fe
	gen b_fe_ci_lo = b_fe - 1.96*se_fe
	gen b_pd_ci_hi = b_pd + 1.96*se_pd
	gen b_pd_ci_lo = b_pd - 1.96*se_pd

	***create ceiling and floor for SEs (for sake of data visualization)
	replace b_fe_ci_hi = .69 if b_fe_ci_hi>.69
	replace b_pd_ci_hi = .69 if b_pd_ci_hi>.69
	replace b_fe_ci_lo = -.09 if b_fe_ci_lo<-.09
	replace b_pd_ci_lo = -.09 if b_pd_ci_lo<-.09	

	***customize clock positions so marker labels do not overlap 
	gen pos = 3
	replace pos = 2 if inlist(pheno, "asthma")
	replace pos = 4 if inlist(pheno, "ever_smk", "phys_act", "birth")
	replace pos = 9 if inlist(pheno, "drinks", "sat_fam", "copd", "health", "open", "sat_fin") // ,
	replace pos = 10 if inlist(pheno, "dep", "aer", "open", "sat_fam")
	
	/*
	replace pos = 5 if inlist(pheno,"health")
	replace pos = 6 if inlist(pheno, "cig_day")
	replace pos = 10 if inlist(pheno, "alch", "asthma", "swb", "dust", "relig")
	*/
	
	*replace pos = 9 if inlist(pheno, "aer", "neuro")

	***estimate correlation and regression coefficient 
	corr b_pd b_fe
	local rho = round(`r(rho)', .001)		
	local rho : display %04.3f `rho'
	
	regress b_pd b_fe
	local slope = round(_b[b_fe], .001)
	local slope : display %03.2f `slope'
	local lo = -.1/`slope'
	local hi = .7/`slope'

	format b_pd b_fe %02.1f
	
	***create scatterplot
	twoway scatter b_pd b_fe, ///
					mlcolor(ebblue%63) ///
					mlwidth(medium) ///   
					mfcolor(ebblue%32) ///
					msize(medlarge) ///
					xline(0, lcolor(black) lpattern(solid)) ///
					yline(0, lcolor(black) lpattern(solid)) ///
					legend(off) ///
					xscale(range(-.1 .7) lstyle(none)) ///
					yscale(range(-.1 .7) lstyle(none)) ///
					title(" " "{bf:B.}" "50% Two Genotype" "Sample", position(12) size(large)) ///	
					subtitle("[1 Rep]" " ", position(12) size(medium)) ///			
					xtitle(" ") ///				
					graphregion(margin(-1 -1 -1 -1)) ///
					plotregion(margin(-1 3 -1 -1)) ///
					ylabel(-.1(.2).7, labcolor(white) labsize(medlarge) notick) ///
					xlabel(-.1(.2).7, labsize(medlarge) notick) ///	
					|| ///
		   scatter b_pd b_fe, ///
					mfcolor(none) mlcolor(none) ///
					mlabel(pheno) mlabcolor(ebblue) mlabvposition(pos) mlabsize(small) mlabangle(315) mlabgap(tiny) ///
					|| ///
		   function y=x, ///
					range(-.1 .7) lpattern(dash) lcolor(black) ///
					|| ///
		   function y = `slope' * x, ///
					range(`lo' `hi') lcolor(maroon%50) ///
					|| ///					
		   rcap b_fe_ci_lo b_fe_ci_hi b_pd if pheno=="bmi", horizontal color(ebblue%32) ///
					|| ///
		   rcap b_pd_ci_lo b_pd_ci_hi b_fe if pheno=="bmi", color(ebblue%32) ///											
		   text(-.1 .5 "{it:r} = `rho'" "β = `slope'", size(medlarge)) ///
		   saving("${figure}/temp/scatter_drop_sibling.gph", replace)

	keep pheno b_pd se_pd
	rename *pd *pd_PANELB		   
	tempfile hold1
	save `hold1', replace

restore

***********************************************************************************
*** NON-OVERLAPPING SAMPLES
***********************************************************************************

matrix drop _all

preserve 

	global rowname ""

	foreach var of global stub {

		global rowname "$rowname `var'"

		quietly reghdfe out_`var' pgi_`var' ///
						if dna_sib==1 ///
						, absorb(idpub)
		local b_fe = _b[pgi_`var']
		local se_fe = _se[pgi_`var']
		local m_2g =  `e(N)'/2
		
		quietly sum rho_pgi_`var'
		quietly phendiff out_`var' pgi_`var' ///
				 if dna_sib==0 ///
				 , fam(idpub) ///
				 rho(`r(mean)') ///
				 delta(2104)
		
		mat betas = nullmat(betas) \ `b_fe', `se_fe', `r(beta)', `r(se)', `m_2g', `r(n)'		
	}			

	mat rownames betas = $rowname			
	mat colnames betas = "b_fe" "se_fe" "b_pd" "se_pd" "m_2g" "n_1g"

	svmat2 betas, names(col) rnames(pheno) 

	keep b_fe se_fe b_pd se_pd m_2g n_1g pheno 
	keep if pheno!=""
	
	gen b_fe_ci_hi = b_fe + 1.96*se_fe
	gen b_fe_ci_lo = b_fe - 1.96*se_fe

	gen b_pd_ci_hi = b_pd + 1.96*se_pd
	gen b_pd_ci_lo = b_pd - 1.96*se_pd

	replace b_fe_ci_hi = .69 if b_fe_ci_hi>.69
	replace b_pd_ci_hi = .69 if b_pd_ci_hi>.69

	replace b_fe_ci_lo = -.09 if b_fe_ci_lo<-.09
	replace b_pd_ci_lo = -.09 if b_pd_ci_lo<-.09	
	
	gen pos = 3
	replace pos = 4 if inlist(pheno, "birth", "swb", "open", "health")
	replace pos = 9 if inlist(pheno, "sat_fin", "sat_fam", "drinks", "phys_act", "cig_day", "copd", "asthma")
	replace pos = 9 if inlist(pheno, "dep")	
	replace pos = 10 if inlist(pheno, "relig")
	
	corr b_pd b_fe
	local rho = round(`r(rho)', .001)	
	local rho : display %04.3f `rho'

	regress b_pd b_fe
	local slope = round(_b[b_fe], .001)
	local slope : display %03.2f `slope'
	local lo = -.1/`slope'
	local hi = .7/`slope'

	format b_pd b_fe %02.1f	
	
	***create scatterplot
	twoway scatter b_pd b_fe, ///
					mlcolor(ebblue%63) ///
					mlwidth(medium) ///   
					mfcolor(ebblue%32) ///
					msize(medlarge) ///
					xline(0, lcolor(black) lpattern(solid)) ///
					yline(0, lcolor(black) lpattern(solid)) ///
					legend(off) ///
					xscale(range(-.1 .7) lstyle(none)) ///
					yscale(range(-.1 .7) lstyle(none)) ///
					title(" " "{bf:C.}" "One Genotype" "Sample", position(12) size(large)) ///		
					subtitle(" " " ", position(12) size(medium)) ///								
					xtitle(" ") ///			
					graphregion(margin(-1 -1 -1 -1)) ///
					plotregion(margin(-1 -1 -1 -1)) ///
					ylabel(-.1(.2).7, labcolor(white) labsize(medlarge) notick) ///
					xlabel(-.1(.2).7, labsize(medlarge) notick) ///								
					|| ///
		   scatter b_pd b_fe, ///
					mfcolor(none) mlcolor(none) ///
					mlabel(pheno) mlabcolor(ebblue) mlabvposition(pos) mlabsize(small) mlabangle(315) mlabgap(tiny) ///
					|| ///
		   function y = x, ///
					range(-.1 .7) lpattern(dash) lcolor(black) ///
					|| ///
		   function y = `slope' * x, ///
					range(`lo' `hi') lcolor(maroon%50) ///		
					|| ///					
		   rcap b_fe_ci_lo b_fe_ci_hi b_pd if pheno=="bmi", horizontal color(ebblue%32) ///
					|| ///
		   rcap b_pd_ci_lo b_pd_ci_hi b_fe if pheno=="bmi", color(ebblue%32) ///						
		   text(-.1 .5 "{it:r} = `rho'" "β = `slope'", size(medlarge)) ///
		   saving("${figure}/temp/scatter_no_overlap.gph", replace)			
		   
	keep pheno b_pd se_pd m_2g n_1g
	rename *pd *pd_PANELC		   
	tempfile hold2
	save `hold2', replace
		   
restore  

***********************************************************************************
*** ONE RANDOM SIBLING SAMPLE (1000 REPS)
***********************************************************************************

matrix drop _all
	
global rowname ""

foreach var of global stub {

	display "`var'"
	
	global rowname "$rowname `var'"

	quietly reghdfe out_`var' pgi_`var' ///
					if dna_sib==1, absorb(idpub)
	local b_fe = _b[pgi_`var']
	local se_fe = _se[pgi_`var']
		
	forval i = 1/${fig1C_reps} {
		preserve
			quietly {
				***randomly select a focal sib
				gen rand = rnormal() if dna_sib==1
				egen fam_mn_rand = mean(rand), by(idpub)
				gen focal = rand>fam_mn_rand
				drop rand fam_mn_rand
				
				***
				replace pgi_`var'=. if focal==0		
				quietly sum rho_pgi_`var'
				quietly phendiff out_`var' pgi_`var' ///
						 if dna_sib==1 ///
						 , fam(idpub) ///
						 rho(`r(mean)') ///
						 delta(2104)
				local se2 = (`r(se)')^2				
				matrix `var' = nullmat(`var') \ `r(beta)', `se2' 
				
			}
		restore	
	}

	***average sibling correlations within phenotype
	capture matrix drop row col avg
	
	mat row = J(rowsof(`var'),1,1)
	mat col = row'*`var'
	mat avg=col/rowsof(`var')
	
	*matrix list `var'
	*matrix list avg
	
	local b_pd = avg[1,1]
	local se_pd = (avg[1,2])^.5
	
	mat betas = nullmat(betas) \ `b_fe', `se_fe', `b_pd', `se_pd'
}			

mat rownames betas = $rowname			
mat colnames betas = "b_fe" "se_fe" "b_pd" "se_pd"
	
***
svmat2 betas, names(col) rnames(pheno) 

keep b_fe se_fe b_pd se_pd pheno
keep if pheno!=""

gen b_fe_ci_hi = b_fe + 1.96*se_fe
gen b_fe_ci_lo = b_fe - 1.96*se_fe

gen b_pd_ci_hi = b_pd + 1.96*se_pd
gen b_pd_ci_lo = b_pd - 1.96*se_pd

replace b_fe_ci_hi = .69 if b_fe_ci_hi>.69
replace b_pd_ci_hi = .69 if b_pd_ci_hi>.69

replace b_fe_ci_lo = -.09 if b_fe_ci_lo<-.09
replace b_pd_ci_lo = -.09 if b_pd_ci_lo<-.09	

gen pos = 3
replace pos = 4 if inlist(pheno, "birth", "swb", "sat_fam")
replace pos = 9 if inlist(pheno, "relig", "copd", "cig_day")
replace pos = 9 if inlist(pheno, "extra", "cog", "asthma", "open", "sat_fin")
replace pos = 9 if inlist(pheno, "ever_smk", "lonely")
replace pos = 10 if inlist(pheno, "phys_act", "extra", "aer", "drinks")

corr b_pd b_fe
local rho = round(`r(rho)', .001)
local rho : display %04.3f `rho'

regress b_pd b_fe
local slope = round(_b[b_fe], .001)
local slope : display %03.2f `slope'
local lo = -.1/`slope'
local hi = .7/`slope'

format b_pd b_fe %02.1f

***create scatterplot
twoway scatter b_pd b_fe, ///
				mlcolor(ebblue%63) ///
				mlwidth(medium) ///   
				mfcolor(ebblue%32) ///
				msize(medlarge) ///			
				xline(0, lcolor(black) lpattern(solid)) ///
				yline(0, lcolor(black) lpattern(solid)) ///
				legend(off) ///
				xscale(range(-.1 .7) lstyle(none)) ///
				yscale(range(-.1 .7) lstyle(none)) ///
				title(" " "{bf:A.}" "50% Two Genotype" "Sample", position(12) size(large)) ///			
				subtitle("[1000 Reps]" " ", position(12) size(medium)) ///							
				xtitle(" ") ///				
				graphregion(margin(-1 -1 -1 -1)) ///
				plotregion(margin(-1 3 -1 -1)) ///
				ylabel(-.1(.2).7, labsize(medlarge) notick) ///							
				xlabel(-.1(.2).7, labsize(medlarge) notick) ///				
				|| ///
	   scatter b_pd b_fe, ///
				mfcolor(none) mlcolor(none) ///
				mlabel(pheno) mlabcolor(ebblue) mlabvposition(pos) mlabsize(small) mlabangle(315) mlabgap(tiny) ///
				|| ///
	   function y=x, ///
				range(-.1 .7) lpattern(dash) lcolor(black) ///
				|| ///
	   function y = `slope' * x, ///
				range(-.1 .7) lcolor(maroon%50) ///
				|| ///					
	   rcap b_fe_ci_lo b_fe_ci_hi b_pd if pheno=="bmi", horizontal color(ebblue%32) ///
				|| ///
	   rcap b_pd_ci_lo b_pd_ci_hi b_fe if pheno=="bmi", color(ebblue%32) ///				
	   text(-.1 .5 "{it:r} = `rho'" "β = `slope'", size(medlarge)) ///			
	   saving("${figure}/temp/scatter_drop_sibling1000.gph", replace) 	   
	   	   
***********************************************************************************
*** COMBINE SCATTERS 
***********************************************************************************	   

graph combine "${figure}/temp/scatter_drop_sibling1000.gph" ///	
			  "${figure}/temp/scatter_drop_sibling.gph" ///	
 			  "${figure}/temp/scatter_no_overlap.gph", ///	
			  ycommon ///
			  col(3) ///
			  l1("Phenotype Differences Estimate             ", size(medlarge)  margin(0 0 0 0) bmargin(0 0 0 0) ) /// margin(zero) bmargin(zero) 
			  b1("Fixed Effects Estimate" " ", margin(zero) bmargin(zero) height(1.2cm)  size(medlarge)) ///
			  imargin(zero)

graph export "${figure}/fig1_${date}.tif", replace width(3300) height(2700)

***********************************************************************************
*** EXPORT SCATTER DATA FOR SI
***********************************************************************************	   

keep pheno b_fe se_fe b_pd se_pd 
rename *pd *pd_PANELA

merge 1:1 pheno using `hold1', keep(3) nogenerate
merge 1:1 pheno using `hold2', keep(3) nogenerate

***label columns
label var pheno "Phenotype"
label var b_fe "Fixed Effects β Estimate"
label var se_fe "Fixed Effects Standard Error"
label var b_pd_PANELA "Phenotype Differences β Estimate (Panel A)"
label var se_pd_PANELA "Phenotype Differences Standard Error (Panel A)"
label var b_pd_PANELB "Phenotype Differences β Estimate (Panel B)"
label var se_pd_PANELB "Phenotype Differences Standard Error (Panel B)"
label var b_pd_PANELC "Phenotype Differences β Estimate (Panel C)"
label var se_pd_PANELC "Phenotype Differences Standard Error (Panel C)"
label var m_2g "Two Genotypes Sample Size (M Sibling Pairs)"
label var n_1g "One Genotype Sample Size (N Sibling Pairs)"

***label rows
replace pheno = "Adventurousness" if pheno=="adv"
replace pheno = "Age at First Birth" if pheno=="birth"
replace pheno = "Age at First Menses" if pheno=="menses"
replace pheno = "Age Voice Deepened" if pheno=="deep"
replace pheno = "Attention Deficit Hyperactivity Disorder" if pheno=="adhd"
replace pheno = "Allergy: Cat" if pheno=="cat"
replace pheno = "Allergy: Dust" if pheno=="dust"
replace pheno = "Allergy: Pollen" if pheno=="pollen"
replace pheno = "Asthma/Eczema/Rhinitis" if pheno=="aer"
replace pheno = "Asthma" if pheno=="asthma" 
replace pheno = "Alcohol Misuse" if pheno=="alcoh"
replace pheno = "Body Mass Index" if pheno=="bmi"
replace pheno = "Cannabis Use" if pheno=="canna"
replace pheno = "Cognitive Empathy" if pheno=="cog_emp"
replace pheno = "Childhood Reading" if pheno=="read" 
replace pheno = "Chronic Obstructive Pulmonary Disease" if pheno=="copd"
replace pheno = "Cigarettes per Day" if pheno=="cig_day"
replace pheno = "Cognitive Performance" if pheno=="cog"
replace pheno = "Delay Discounting" if pheno=="dly_disc"
replace pheno = "Depressive Symptoms" if pheno=="dep"
replace pheno = "Drinks per Week" if pheno=="drinks"
replace pheno = "Educational Attainment" if pheno=="edu"
replace pheno = "Ever Smoker" if pheno=="ever_smk" 
replace pheno = "Extraversion" if pheno=="extra" 
replace pheno = "Life Satisfaction: Family" if pheno=="sat_fam"
replace pheno = "Life Satisfaction: Finance" if pheno=="sat_fin"
replace pheno = "Life Satisfaction: Friend" if pheno=="sat_frnd" 
replace pheno = "Life Satisfaction: Work" if pheno=="sat_job"
replace pheno = "Hayfever" if pheno=="hay" 
replace pheno = "Height" if pheno=="hgt" 
replace pheno = "Highest Math" if pheno=="high_math" 
replace pheno = "Left Out of Social Activity" if pheno=="leftout"
replace pheno = "Loneliness" if pheno=="lonely"
replace pheno = "Morning Person" if pheno=="chrono"
replace pheno = "Migraine" if pheno=="migrn" 
replace pheno = "Narcissism" if pheno=="narci"
replace pheno = "Nearsightedness" if pheno=="near_sgt" 
replace pheno = "Number Ever Born: Men" if pheno=="neb_male" 
replace pheno = "Number Ever Born: Women" if pheno=="neb_fem" 
replace pheno = "Neuroticism" if pheno=="neuro" 
replace pheno = "Openness" if pheno=="open" 
replace pheno = "Physical Activity" if pheno=="phys_act"
replace pheno = "Religious Attendance" if pheno=="relig"
replace pheno = "Risk Tolerance" if pheno=="risk"
replace pheno = "Self-Rated Health" if pheno=="health"
replace pheno = "Self-Rated Math Ability" if pheno=="self_math"
replace pheno = "Subjective Well-Being" if pheno=="swb"

order pheno
sort pheno

export excel using "${table}/tabS4_${date}.xlsx", firstrow(varlabels) replace keepcellfmt
