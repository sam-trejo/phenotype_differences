
do "${syntax}/programs/bky2006.do"
	
***********************************************************************************
*** PANEL A 
***********************************************************************************

use "${full}", clear

collapse (mean) dna, by(died)
keep if inrange(died, 1980, 2020)

format dna %02.1f

twoway connect dna died, ///
	   mlcolor(ebblue*.63) ///
	   mlwidth(medium) ///   
	   mfcolor(ebblue*.63) ///
	   msize(medlarge) ///
	   msymbol(O) ///
	   msize(medium) ///
   	   lwidth(medthick) ///   	   	   
	   lcolor(ebblue*.63) ///   
	   xline(2006.5, lcolor(black)) ///
	   xlabel(1985(10)2020, labsize(medlarge)) ///
	   ylabel(, labsize(medlarge)) ///		   
	   xtick(1980(5)2020) ///
	   title(" {bf:A.} Mortality Selection into Genotyping", position(12) size(huge)) ///	 justification(left)
	   xtitle("Death Year", size(vlarge)) ///   
	   ytitle("Fraction Genotyped", size(vlarge)) ///  orientation(horizontal) 
       legend(order(1 "Graduates" 2 "Siblings") cols(1) position(9) ring(0) size(large)) /// 
	   graphregion(margin(0 17 0 0)) ///
	   saving("${figure}/temp/genotyping.gph", replace) 		

***********************************************************************************
*** PANEL B
***********************************************************************************

use "${analytic}", clear

***histogram
twoway hist diff_imp_span if dna_sib==0 & is_grad==1 & sib_died==1,  ///
			freq color(ebblue%50) lcolor(ebblue) lwidth(medium) ///
			width(2) start(-52) ///
	   || ///
	   hist diff_imp_span if dna_sib==1 & is_grad==1 & sib_died==1, ///
			freq color(maroon%50) lcolor(maroon) lwidth(medium) ///
			width(2) start(-52) /// 
			xlabel(-50(25)50, labsize(medlarge)) ///
			ylabel(, labsize(medlarge)) ///	
			xtick(-50(12.5)50) ///			
			yscale(range(0 170)) ///
			title("{bf:B.} Within-Family Variation in Lifespan  ", position(12) size(huge)) ///			
			xtitle("Graduate Lifespan – Sibling Lifespan", size(vlarge)) ///
			ytitle("{it:N} Sibling Pairs", size(vlarge)) /// orientation(horizontal)
			legend(order(1 "1x Genotype Sample" 2 "2x Genotype Sample") position(2) ring(0) size(medlarge) col(1) bmargin(0 0 0 0)) ///    rowgap(*.15)
			graphregion(margin(0 17 0 0)) ///
			saving("${figure}/temp/hist_span_diff.gph", replace) 

***********************************************************************************
*** PANEL C
***********************************************************************************

*** LOAD WLS DATA
use "${analytic}", clear

*** DROP EXTRA PGI
drop *pgi_menses *pgi_neb_male *pgi_neb_fem *pgi_deep ///
	 *pgi_adhd *pgi_pollen *pgi_risk 
	 
*** REGRESS LIFESPAN ON PGIS
matrix drop _all
global rowname ""

***loop over PGIs
foreach var of varlist pgi* {
		
	global rowname "$rowname `var'"
	
	***run fixed effects regression
	quietly reghdfe imp_span `var' ///
					female c_born c_born2 female_c_born female_c_born2 ///
					if dna_sib==1, absorb(idpub) 
	local b_fe = _b[`var']
	local se_fe = _se[`var']
	local p_fe = 2*ttail(e(df_r),abs(`b_fe'/`se_fe')) 
	local n_fe = round(`e(N)'*.5,.001)
	
	***run phenotype differences regression
	quietly sum rho_`var'
	local rho = `r(mean)'	
	quietly reg res_diff_imp_span x5_`var' if dna_sib==0
	local b_pd = _b[x5_`var']
	local adj = `b_pd'^2 * (1+`rho') / (1-`rho') / (2105)
	local se_pd = (_se[x5_`var']^2 + `adj')^.5	
	local p_pd = 2*ttail(e(df_r),abs(`b_pd'/`se_pd')) 
	local n_pd = round(`e(N)',.001)

	***output regression coefficients to a matrix
	mat betas = nullmat(betas) \ `b_fe', `se_fe', `p_fe', `n_fe', `b_pd', `se_pd', `p_pd', `n_pd'
}			

***name rows and columns of regression coefficient matrix
mat rownames betas = $rowname			
mat colnames betas = "b_fe" "se_fe" "p_fe" "n_fe" "b_pd" "se_pd" "p_pd" "n_pd"

***tranfer regression coefficient matrix to a variable
svmat2 betas, names(col) rnames(pheno) 

***
keep b_fe se_fe p_fe n_fe b_pd se_pd p_pd n_pd pheno
keep if pheno!=""

***generate inverse variance weighted meta-anaylsis of FE and PD estimates
gen var_fe = se_fe^2
gen var_pd = se_pd^2
gen b_meta_iv = ((b_fe/var_fe)+(b_pd/var_pd))/(1/var_fe+1/var_pd)
gen se_meta_iv = (1/(1/var_fe+1/var_pd))^.5
gen p_meta_iv = 2*ttail(n_fe+n_pd,abs(b_meta_iv/se_meta_iv)) 

***
keep pheno p_meta_iv b_meta_iv se_meta_iv
rename *_meta_iv* **
rename p pval
sort pval
drop if pval==.

tempfile hold
save `hold', replace
keep if pheno=="pgi_meta"
tempfile meta
save `meta', replace
use `hold', clear

drop if pheno=="pgi_meta"
bky2006	

append using `meta'
rename qval qval2
bky2006	

order pheno b se pval qval qval2
sort pval

***merge on measurement error multipliers
merge 1:1 pheno using "${data}/me_mult_long_2024_01_05.dta", nogenerate keep(1 3)

***generate disattenuated betas and ses
gen b_mec = b * me_mult if qval<.1
keep pheno b se pval qval qval2 b_mec
order pheno b se pval qval2 b_mec
sort b

gen high = b + 1.96*se
gen low = b - 1.96*se
gen n = _n

***Label Phenotypes
replace pheno = "Adventurousness" if pheno=="pgi_adv"
replace pheno = "Age at First Birth" if pheno=="pgi_birth"
replace pheno = "Allergy: Cat" if pheno=="pgi_cat"
replace pheno = "Allergy: Dust" if pheno=="pgi_dust"
replace pheno = "Asthma/Eczema/Rhini." if pheno=="pgi_aer"
replace pheno = "Asthma" if pheno=="pgi_asthma" 
replace pheno = "Alcohol Misuse" if pheno=="pgi_alcoh"
replace pheno = "{bf:Body Mass Index}" if pheno=="pgi_bmi"
replace pheno = "Cannabis Use" if pheno=="pgi_canna"
replace pheno = "Cognitive Empathy" if pheno=="pgi_cog_emp"
replace pheno = "Childhood Reading" if pheno=="pgi_read" 
replace pheno = "{bf:Chronic Obstructive P.}" if pheno=="pgi_copd"
replace pheno = "Cigarettes per Day" if pheno=="pgi_cig_day"
replace pheno = "Cognitive Performance" if pheno=="pgi_cog"
replace pheno = "Delay Discounting" if pheno=="pgi_dly_disc"
replace pheno = "Depressive Symptoms" if pheno=="pgi_dep"
replace pheno = "Drinks per Week" if pheno=="pgi_drinks"
replace pheno = "Educational Attainment" if pheno=="pgi_edu"
replace pheno = "Ever Smoker" if pheno=="pgi_ever_smk" 
replace pheno = "Extraversion" if pheno=="pgi_extra" 
replace pheno = "Life Satisfaction: Family" if pheno=="pgi_sat_fam"
replace pheno = "Life Satisfaction: Finance" if pheno=="pgi_sat_fin"
replace pheno = "Life Satisfaction: Friend" if pheno=="pgi_sat_frnd" 
replace pheno = "Life Satisfaction: Work" if pheno=="pgi_sat_job"
replace pheno = "Hayfever" if pheno=="pgi_hay" 
replace pheno = "Height" if pheno=="pgi_hgt" 
replace pheno = "Highest Math" if pheno=="pgi_high_math" 
replace pheno = "Left Out of Social A." if pheno=="pgi_leftout"
replace pheno = "Loneliness" if pheno=="pgi_lonely"
replace pheno = "Morning Person" if pheno=="pgi_chrono"
replace pheno = "Migraine" if pheno=="pgi_migrn" 
replace pheno = "Narcissism" if pheno=="pgi_narci"
replace pheno = "Nearsightedness" if pheno=="pgi_near_sgt" 
replace pheno = "Neuroticism" if pheno=="pgi_neuro" 
replace pheno = "Openness" if pheno=="pgi_open" 
replace pheno = "Physical Activity" if pheno=="pgi_phys_act"
replace pheno = "Religious Attend." if pheno=="pgi_relig"
replace pheno = "{bf:Self-Rated Health}" if pheno=="pgi_health"
replace pheno = "Self-Rated Math Abil." if pheno=="pgi_self_math"
replace pheno = "Subjective Well-Being" if pheno=="pgi_swb"
replace pheno = "{bf:Meta}" if pheno=="pgi_meta"

format high low b %02.1f

gen marker = -.499
replace marker = .102 if n<19

twoway rcap high low n if qval<.1, ///
	   horizontal ///
	   lcolor(ebblue*.63) ///
   	   lwidth(medthick) ///   	   	   
	   || ///
	   rcap high low n if qval>=.1, ///
	   horizontal ///
	   lcolor(ebblue*.32) ///
   	   lwidth(medthick) ///   	   	   
	   || ///
	   scatter n b if qval<.1, ///
	   mlcolor(ebblue*.63) ///
	   mlwidth(medium) ///   
	   mfcolor(ebblue*.63) ///
	   msize(large) ///
	   msymbol(o) ///	   
	   || ///   
	   scatter n b if qval>=.1, ///
	   mcolor(ebblue*.32) ///
	   mlwidth(medium) ///   
	   msize(large) ///
	   msymbol(o) ///	   
	   || /// 
	   scatter n marker, ///
		msymbol(none) ///
		mlabel(pheno) mlabcolor(black) mlabsize(medsmall) ///		
		xline(0) ///
		xtick(-.5(.1).5) ///
	    yscale(lstyle(none)) ///
		xlabel(-.4(.2).4, labsize(medsmall)) ///
		ylabel(, nolabels noticks) ///
		xtitle("Lifespan (Years)", size(medlarge)) ///
		ytitle("") ///
		legend(off) ///
		title("{bf:C.} Effects of Polygenic Scores on Lifespan", position(12) size(large)) ///		
		graphregion(margin(0 0 2 2)) ///
		saving("${figure}/temp/lifespan.gph", replace)
		
***********************************************************************************
*** COMBINE PANELS A AND B, and then add PANEL C on the right
***********************************************************************************

graph combine ///
	  "${figure}/temp/genotyping.gph" ///
	  "${figure}/temp/hist_span_diff.gph", ///
	  cols(1) ///
 	  saving("${figure}/temp/selection.gph", replace)

graph combine ///
	  "${figure}/temp/selection.gph" ///
	  "${figure}/temp/lifespan.gph", ///
	  cols(2) ///
	  graphregion(margin(zero))	
		
graph export "${figure}/fig2_${date}.tif", replace height(2700) width(3300)
		
***********************************************************************************
*** EXPORT BETA ESTIMATES FOR SI
***********************************************************************************

label var pheno "Polygenic Score"
label var b "β Estimate (Outcome: Lifespan)"
label var se "Standard Error"
label var pval "p-Value"
label var qval "q-Value (including meta-PGS)"
label var qval2 "q-Value (excluding meta-PGS)"
label var b_mec "Disattenuated β Estimate "

replace pheno = "Asthma/Eczema/Rhinitis" if pheno=="Asthma/Eczema/Rhini."
replace pheno = "Body Mass Index" if pheno=="{bf:Body Mass Index}"
replace pheno = "Chronic Obstructive Pulmonary Disease" if pheno=="{bf:Chronic Obstructive P.}"
replace pheno = "Left Out of Social Activity" if pheno=="Left Out of Social A."
replace pheno = "Religious Attendance" if pheno=="Religious Attend."
replace pheno = "Self-Rated Health" if pheno=="{bf:Self-Rated Health}"
replace pheno = "Meta" if pheno=="{bf:Meta}"

keep pheno b se pval qval qval2 b_mec
order pheno b se pval qval qval2 b_mec
sort pheno

export excel using "${table}/tabS3_${date}.xlsx", firstrow(varlabels) replace
