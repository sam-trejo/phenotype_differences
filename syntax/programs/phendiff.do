
***********************************************************************************
*** PHENOTYPE DIFFERENCES REGRESSION PROGRAM
***********************************************************************************

capture program drop phendiff

program define phendiff, rclass
syntax varlist(min=2 max=2) [if] [, Fam(varname) RHO(real 0.5) RESidualize(varlist) DELTa(int -1) DESTroy] 

	if "`destroy'"=="" {
		preserve	
	}
			local out: word 1 of `varlist'
		 
			local g: word 2 of `varlist'
				
			quietly {
				
				if "`if'"!="" {
					keep `if'
				}
				
				***drop families without exactly two siblings			
				egen fam_n = count(`fam'), by(`fam')
				keep if fam_n==2 			
				
				***drop families without exactly two phenotypes observed
				egen fam_n_`out' = count(`out'), by(`fam')
				keep if fam_n_`out'==2 
				
				***drop families without exactly one genotype observed
				egen fam_n_`g' = count(`g'), by(`fam')
				keep if fam_n_`g'==1 						
				
				***generate sibling phenotype difference			
				egen fam_mn_`out' = mean(`out'), by(`fam')
				gen diff_`out' = (`out' - fam_mn_`out') * 2
				
				***generate rho-transformed genetic predictor		
				gen `g'_X_rho = (1 - `rho') * `g'

				***run phenotype differences regression (without covariate residualization)
				if "`residualize'"=="" {
					reg diff_`out' `g'_X_rho	
				}	
				
				***run phenotype differences regression (with covariate residualization)
				if "`residualize'"!="" {
					foreach var of varlist `residualize' {
				
						egen fam_n_`var' = count(`var'), by(`fam')
						keep if fam_n_`var'==2 	// ***drop families without exactly two of each covariate
						
						egen fam_mn_`var' = mean(`var'), by(`fam')
						gen diff_`var' = (`var' - fam_mn_`var') * 2 // ***generate sibling covariate difference
						replace `var' =  diff_`var' 			
					}
					
					reg diff_`out' `residualize'
					predict res_diff_`out' if e(sample), residuals
					
					reg res_diff_`out' `g'_X_rh
				}
			}	

			***store results		
			local beta = _b[`g'_X_rho]
			local n = `e(N)'
			local se = _se[`g'_X_rho]
			local se_prefix = "SE   = "		
			
			***adjust standard errors using the delta method
			if `delta'>0 {
				local adj = `beta'^2 * (1+`rho')^2 / (`delta')
				local se = (_se[`g'_X_rho]^2 + `adj')^.5	
				local se_prefix = "SE   ≈ "
			}					
			
			***return results	
			return local n `n'			
			return local se `se'
			return local beta `beta'				
			
			***display results
			local n_display = "Phenotype differences regression model using " + strofreal(`n', "%15.0fc") + " sibling pairs"
			local beta_display = "Beta = " + strofreal(`beta', "%15.3f")
			local se_display = "`se_prefix'" + strofreal(`se', "%15.3f")
			display "—————————————————————————————————————————————————————————————————"	
			display "`n_display'"
			display "—————————————————————————————————————————————————————————————————"	
			display "`beta_display'"
			display "`se_display'"			
end
