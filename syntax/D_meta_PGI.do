
***********************************************************************************
*** 
***********************************************************************************

use "${analytic_temp}", clear

if "${enet}"=="YES" {
	preserve
		drop pgi_menses pgi_neb_male pgi_neb_fem pgi_deep ///
			 pgi_adhd pgi_pollen pgi_risk
			
		drop if dna==0		
		egen mn_born = mean(born), by(idpub)		
		drop if born>mn_born		
		codebook idpub		
			
		quietly elasticnet linear imp_span (female c_born c_born2 female_c_born female_c_born2 pc*) pgi*, ///
						   alpha(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1) ///
						   rseed(19146) crossgrid(union) alllambdas  
		display "`e(alpha_sel)'"						   
	restore

	matselrc e(b) beta, r(1) c(`e(othervars_sel)')
	matewmf beta abs_beta, f(abs)
	mata : st_matrix("tot_abs_beta", rowsum(st_matrix("abs_beta")))
	matrix scale_beta = (1/tot_abs_beta[1,1])*beta

	local n = `e(k_nonzero_cv)' - `e(k_nonzero_serule)' 
	display `n'
	forval i=1/`n' {
		matrix scale_beta[1,`i'] = round(scale_beta[1,`i'], .01)

	}
	matrix list scale_beta
}

***********************************************************************************
*** 
***********************************************************************************

gen pgi_meta =  -.25*pgi_bmi      +  .11*pgi_edu    + -.11*pgi_ever_smk + /// 
				-.12*pgi_hgt      +  .03*pgi_narci  + -.15*pgi_near_sgt + ///
				 .01*pgi_read     + -.01*pgi_extra  +  .09*pgi_sat_fam  + ///
				 .01*pgi_sat_frnd +  .10*pgi_health

***standardize pgi (to distribution of graduates)
foreach var of varlist pgi* {
	sum pgi_meta if rtype=="g"
	replace pgi_meta = (pgi_meta - `r(mean)')/`r(sd)'
}

save "${analytic_temp_plus_meta}", replace
			   	

