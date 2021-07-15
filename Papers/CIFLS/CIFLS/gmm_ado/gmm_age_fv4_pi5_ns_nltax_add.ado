program gmm_age_fv4_pi5_ns_nltax_add
	version 11
	syntax varlist if, at(name)
	quietly {
		local eps_string ""
		forvalues i=1/31 {
			local eps_string "`eps_string' eps_`i'"
			}

		tempvar psi_c_lam psi_c_w1 psi_c_w2 ///
				psi_h1_lam psi_h1_w1 psi_h1_w2 ///
				psi_h2_lam psi_h2_w1 psi_h2_w2 /// 
				psi_h_lam psi_h_w1 psi_h_w2 ///
				k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 ///
				pi denom `eps_string'
			
		matrix debug_est = `at'

		local eta_c_p	= `at'[1,1]
		local eta_c_w1	= `at'[1,2]
		local eta_h1_w1	= `at'[1,3]
		local eta_c_w2	= `at'[1,4]
		local eta_h1_w2	= `at'[1,5]
		local eta_h2_w2	= `at'[1,6]
		local a0		= `at'[1,7]

		local eta_h1_p	= -`eta_c_w1'*${c_over_y}
		local eta_h2_p	= -`eta_c_w2'*${c_over_yw}
		local eta_h2_w1 = `eta_h1_w2'*${y_over_yw}

		local s_omega 	= $s_omega
		
		gen `pi'	= (pi_hat * `a0') `if' 

		gen `psi_c_lam'=.
		gen `psi_c_w1'=.
		gen `psi_c_w2'=.
		gen `psi_h1_lam'=.
		gen `psi_h1_w1'=.
		gen `psi_h1_w2'=.
		gen `psi_h2_lam'=.
		gen `psi_h2_w1'=.
		gen `psi_h2_w2'=.

		noisily	mata: calc_invQR()

		gen `psi_h_lam' = s_hat*`psi_h1_lam' + (1-s_hat)*`psi_h2_lam'
		gen `psi_h_w1' 	= s_hat*`psi_h1_w1' + (1-s_hat)*`psi_h2_w1'
		gen `psi_h_w2' 	= s_hat*`psi_h1_w2' + (1-s_hat)*`psi_h2_w2'
		gen `denom'  	= `psi_c_lam' - (1-tau)*(`pi')*`psi_h_lam'
		
		gen `k1' = `psi_c_w1'
		gen `k2' = `psi_c_w2'
		gen `k3' = `psi_c_w1' + `psi_c_lam'*((1-tau)*(`pi')*(s_hat+`psi_h_w1')-`psi_c_w1')/`denom'
		gen `k4' = `psi_c_w2' + `psi_c_lam'*((1-tau)*(`pi')*((1-s_hat)+`psi_h_w2')-`psi_c_w2')/`denom'
		
		gen `k5' = 1 + `psi_h1_w1'
		gen `k6' = `psi_h1_w2'
		gen `k7' = 1 + `psi_h1_w1' + `psi_h1_lam'*((1-tau)*(`pi')*(s_hat+`psi_h_w1') - `psi_c_w1')/`denom'
		gen `k8' = `psi_h1_w2' + `psi_h1_lam'*((1-tau)*(`pi')*((1-s_hat)+`psi_h_w2') - `psi_c_w2')/`denom'

		gen `k9' = `psi_h2_w1'
		gen `k10' = 1 + `psi_h2_w2'
		gen `k11' = `psi_h2_w1' + `psi_h2_lam'*((1-tau)*(`pi')*(s_hat+`psi_h_w1') - `psi_c_w1')/`denom'
		gen `k12' = 1 + `psi_h2_w2' + `psi_h2_lam'*((1-tau)*(`pi')*((1-s_hat)+`psi_h_w2') - `psi_c_w2')/`denom'

	
		#delimit;
		//Variance of Consumption
		local j=1;
		gen `eps_`j'' 	= duc2 - ((`k1')^2 * m_s_u1*var_dQ //Var(q1)
										+ (`k2')^2 * m_s_u2*var_dQ //var(q2)
										+ (`k3')^2 * m_s_v1*var_dP //Var(p1)
										+ (`k4')^2 * m_s_v2*var_dP //Var(p2)
										+ (`k1'*`k2') * m_r_u1u2*var_dQ //covar(q1,q2)
										+ (`k3'*`k4') * m_r_v1v2*var_dP //covar(p1,p2)   covariance moments very similar
								+ 2*`s_omega'); //what is s_omega, No longer multiplying by 2.  Checked -10/6 
								
		//		gen `eps_`j'' 	= 	duc2 - 2*((`k1')^2*m_s_u1 + (`k2')^2*m_s_u2 
		// 						+ (`k3')^2*m_s_v1 + (`k4')^2*m_s_v2 
		// 						+ `k1'*`k2'*(2*m_r_u1u2)
		// 						+ `k3'*`k4'*(2*m_r_v1v2)
		// 						+ `s_omega');

		//Variance of Male Income
		local j=`j'+1;
		gen `eps_`j'' 	=	duy2 - (  (`k5')^2 * (a *(1/3)*m_s_u1 +  (1-a*(2/3))*m_s_u1  +  a*(1/3)*m_s_u1  +  (1-a*(2/3))*m_s_u1)  //Var(qbar(1,t))
									+ (`k6')^2 * (a *(1/3)*m_s_u2 +  (1-a*(2/3))*m_s_u2  +  a*(1/3)*m_s_u2  +  (1-a*(2/3))*m_s_u2) //Var(qbar(2,t))
									+ (`k7')^2 * ((1/3)*m_s_v1  +  m_s_v1  +  (1/3)*m_s_v1) //Var(pbar(1,t))
									+ (`k8')^2 * ((1/3)*m_s_v2  +  m_s_v2  +  (1/3)*m_s_v2) //Var(pbar(2,t))
									+ (`k5'*`k6') * (a *(1/3)*m_r_u1u2 +  (1-a*(2/3))*m_r_u1u2  +  a*(1/3)*m_r_u1u2  +  (1-a*(2/3))*m_r_u1u2) //Covar(qbar(1,t), qbar(2,t))  //covariance moments are very identical to variance
									+ (`k7'*`k8') * ((1/3)*m_r_v1v2  +  m_r_v1v2  +  (1/3)*m_r_v1v2) //Covar(pbar(1,t), pbar(2,t))  
						+ s_me_y);  // Checked 10/6
						
		// 		gen `eps_`j'' 	=	duy2 - 2*((`k5')^2*m_s_u1 + (`k6')^2*m_s_u2
		// 						+ (`k7')^2*m_s_v1 + (`k8')^2*m_s_v2
		// 						+ `k5'*`k6'*(2*m_r_u1u2)
		// 						+ `k7'*`k8'*(2*m_r_v1v2)
		// 						+ s_me_y);

		//Variance of Female Income		
		local j=`j'+1;
		gen `eps_`j'' 	=	duyw2 -   ((`k9')^2 * (a *(1/3)*m_s_u1 +  (1-a*(2/3))*m_s_u1  +  a*(1/3)*m_s_u1  +  (1-a*(2/3))*m_s_u1) //Var(qbar(1,t))
									+ (`k10')^2 * (a *(1/3)*m_s_u2 +  (1-a*(2/3))*m_s_u2  +  a*(1/3)*m_s_u2  +  (1-a*(2/3))*m_s_u2) //Var(qbar(2,t))
									+ (`k11')^2 * ((1/3)*m_s_v1  +  m_s_v1  +  (1/3)*m_s_v1) //Var(pbar(1,t))
									+ (`k12')^2 * ((1/3)*m_s_v2  +  m_s_v2  +  (1/3)*m_s_v2) //Var(pbar(2,t))
									+ (`k9'*`k10') * (a *(1/3)*m_r_u1u2 +  (1-a*(2/3))*m_r_u1u2  +  a*(1/3)*m_r_u1u2  +  (1-a*(2/3))*m_r_u1u2) //Cov(qbar(1,t), qbar(2,t))   //covariance moments are very identical to variance
									+ (`k11'*`k12') * ((1/3)*m_r_v1v2  +  m_r_v1v2  +  (1/3)*m_r_v1v2)  //Cov(pbar(1,t), pbar(2,t))
						+ s_me_yw_cons); //Checked 10-6
						
		// 		gen `eps_`j'' 	=	duyw2 - ((`k9')^2*(m_s_u1)*(2-(r_u1)^2*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
		// 						+ (`k10')^2*(m_s_u2)*(2-(r_u2)^2*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
		// 						+ (`k11')^2*(m_s_v1)*(2-(r_v1)^2*z_eta*inv_mills) 
		// 						+ (`k12')^2*(m_s_v2)*(2-(r_v2)^2*z_eta*inv_mills)
		// 						+ 2*(`k9')*(`k10')*(2*m_r_u1u2 - r_u1*r_u2*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
		// 						+ 2*(`k11')*(`k12')*(2*m_r_v1v2 - r_v1*r_v2*z_eta*inv_mills)
		// 						+ (2)*s_me_yw_cons);

		//Covariance between Male income and lagged male income. 
		local j=`j'+1;
		gen `eps_`j'' 	=	duy_lag - ( `k5'^2 * ((-1/3) * a * l2_s_u1  -  (1-(2/3)*a)*l2_s_u1)
									+   `k6'^2 * ((-1/3) * a * l2_s_u2  -  (1-(2/3)*a)*l2_s_u2)
									+  2*`k5'*`k6' * ((-1/3) * a * l2_r_u1u2  -  (1-(2/3)*a)*m_s_u2)
									- s_me_y) ;
									
	// 		gen `eps_`j'' 	=	duy_lag - (-(`k5')^2*(l2_s_u1) -(`k6')^2*(l2_s_u2)
	// 						- 2*(`k5')*(`k6')*l2_r_u1u2
	// 						- s_me_y);

		//Covariance between Female Income and Lagged female income. 
		local j=`j'+1;
		gen `eps_`j'' 	=	duyw_lag - 	( (`k9')^2  *  (a*(-1/3)*l2_s_u1  -  (1-a*(2/3))*l2_s_u1)  //Cov(qbar(1,t), qbar(1,t-2))  
										+ (`k10')^2  *  (a*(-1/3)*l2_s_u2  -  (1-a*(2/3))*l2_s_u2)  //Cov(qbar(2,t), qbar(2,t-2))
										+ 2*(`k9'*`k10')  *  (a*(-1/3)*l2_r_u1u2  -  (1-a*(2/3))*l2_r_u1u2)
										- s_me_yw_cons); //Checked 10-6
										
		// 		gen `eps_`j'' 	=	duyw_lag - (-(`k9')^2*(l2_s_u1)*(1-(r_u1)^2*l2_z_eta*l2_inv_mills)
		// 						-(`k10')^2*(l2_s_u2)*(1-(r_u2)^2*l2_z_eta*l2_inv_mills)
		// 						-2*(`k9')*(`k10')*(l2_r_u1u2-r_u1*r_u2*l2_z_eta*l2_inv_mills)
		// 						- s_me_yw_cons);

		//Covariance between male wage and consumption. 
		local j=`j'+1;
		gen `eps_`j''	=	duw_duc - ( `k1' * m_s_u1 * cov_dQ_bardQ   //Cov(qbar(1,t), q(1,t))
									+   `k2' * m_r_u1u2 * cov_dQ_bardQ //Cov(qbar(1,t), q(2,t))
									+   `k3' * m_s_v1 * cov_dP_bardP   //Cov(pbar(1,t), p(1,t))
									+   `k4' * m_r_v1v2 * cov_dP_bardP); //Cov(pbar(1,t), p(2,t))  //Checked 10-6
									
		// 		gen `eps_`j''	=	duw_duc - 2*((`k1')*(m_s_u1)
		// 						+ `k2'*(m_r_u1u2) 
		// 						+ `k3'*(m_s_v1)
		// 						+ `k4'*(m_r_v1v2));


		//Covariance between Male Wage and Male Income 
		local j=`j'+1;
		gen `eps_`j''	= 	duw_duy - (`k5' * m_s_u1 * var_bardQ    //Var(qbar(1,t))
									+  `k6' * m_r_u1u2 * var_bardQ  //Cov(qbar(1,t), qbar(2,t))
									+  `k7' * m_s_u1 * var_bardP    //Var(pbar(1,t))
									+  `k8' * m_r_v1v2 * var_bardP  //Cov(pbar(1,t), pbar(2,t))
									+  s_me_y - cov_me_y_h); //Checked 10-6
									
		// 		gen `eps_`j''	= 	duw_duy - 2*(`k5'*(m_s_u1) + `k6'*(m_r_u1u2)
		// 						+ `k7'*(m_s_v1) 
		// 						+ `k8'*(m_r_v1v2)
		// 						+ s_me_y - cov_me_y_h);
		

		//Covariance between Male Wage and Female Income 
		local j=`j'+1;
		gen `eps_`j''	=	duw_duyw - (`k9' * m_s_u1 * var_bardQ    //Var(qbar(1,t))
									 +  `k10' * m_r_u1u2 * var_bardQ //Cov(qbar(1,t), qbar(2,t))
									 +  `k11' * m_s_v1 * var_bardP   //Var(pbar(1,t))
									 +  `k12' * m_r_v1v2 * var_bardP); //Cov(pbar(1,t), pbar(2,t))  //no measurement error? How come?	Checked 10-6
		
		//gen `eps_`j''	=	duw_duyw - ((`k9')*(2*m_s_u1 - (r_u1)^2*z_eta*inv_mills - (l2_r_u1)^2*l2_z_eta*l2_inv_mills)
						//+ `k10'*(2*m_r_u1u2 - (r_u1*r_u2)*z_eta*inv_mills - (l2_r_u1*l2_r_u2)*l2_z_eta*l2_inv_mills)
						//+ `k11'*(2*m_s_v1-(r_v1)^2*z_eta*inv_mills)
						//+ `k12'*(2*m_r_v1v2 - (r_v1*r_v2)*z_eta*inv_mills));	
		
		//Covariance between Female Wage and Consumption
		local j=`j'+1;
		gen `eps_`j''	=	duww_duc - (`k1' * m_r_u1u2 * cov_dQ_bardQ    //Cov(qbar(2,t), q(1,t))
									+   `k2' * m_s_u2 * cov_dQ_bardQ     //Cov(qbar(2,t), q(2,t))
									+   `k3' * m_r_v1v2 * cov_dP_bardP //Cov(pbar(2,t), p(1,t))
									+   `k4' * m_s_v2 * cov_dP_bardP );   //Cov(pbar(2,t), p(2,t))  Checked 10-6
		
		//gen `eps_`j''	=	duww_duc - ((`k1')*(2*m_r_u1u2 - r_u1*r_u2*z_eta*inv_mills - l2_r_u1*l2_r_u2*l2_z_eta*l2_inv_mills)
						//+ `k2'*(2*m_s_u2 - (r_u2)^2*z_eta*inv_mills - (l2_r_u2)^2*l2_z_eta*l2_inv_mills)
						//+ `k3'*(2*m_r_v1v2 - (r_v1*r_v2)*z_eta*inv_mills)
						//+ `k4'*(2*m_s_v2 - (r_v2)^2*z_eta*inv_mills));
		
		//Covariance between Female Wage and Male Income 
		local j=`j'+1;
		gen `eps_`j''	=	duww_duy - (`k5' * m_r_u1u2 * var_bardQ   //Cov(qbar(2,t), qbar(1,t))
									+   `k6' * m_s_u2 * var_bardQ   //Var(qbar(2,t))
									+   `k7' * m_r_v1v2 * var_bardP //Cov(pbar(2,t), pbar(1,t))
									+   `k8' * m_s_v2 * var_bardP ); //Var(pbar(2,t))  Checked 10-6   //no measurement error originally Checked 10-6
		
		//gen `eps_`j''	=	duww_duy - (`k5'*(2*m_r_u1u2 - (r_u1*r_u2)*(z_eta*inv_mills) - (l2_r_u1*l2_r_u2)*(l2_z_eta*l2_inv_mills))
						//+ `k6'*(2*(m_s_u2)-(r_u2)^2*(z_eta*inv_mills) - (l2_r_u2)^2*l2_z_eta*l2_inv_mills)
						//+ `k7'*(2*m_r_v1v2 - r_v1*r_v2*z_eta*inv_mills)
						//+ `k8'*(2*(m_s_v2)-(r_v2)^2*z_eta*inv_mills));
		
		///////////////////////////////////////////
		//Covariance between female wage and female income
		
		local j=`j'+1;
		gen `eps_`j''	=	duww_duyw - (`k9' * m_r_u1u2 * var_bardQ    //Cov(qbar(2,t), qbar(1,t))
									+   `k10' * m_s_u2 * var_bardQ     //Var(qbar(2,t))
									+   `k11' * m_r_v1v2 * var_bardP  //Cov(pbar(2,t), pbar(1,t))
									+   `k12' * m_s_v2 * var_bardP //Var(pbar(2,t))
						+ (2)*s_me_yw_cons - (2)*cov_me_yw_hw_cons); //measurement error included   Checked 10-6
						

	// 		gen `eps_`j''	=	duww_duyw - (`k9'*(2*m_r_u1u2 - (r_u1*r_u2)*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
	// 						+ `k10'*(m_s_u2)*(2-(r_u2)^2*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
	// 						+ `k11'*(2*m_r_v1v2 - (r_v1*r_v2)*z_eta*inv_mills)
	// 						+ `k12'*(m_s_v2)*(2-(r_v2)^2*z_eta*inv_mills) 
	//  						+ (2)*s_me_yw_cons - (2)*cov_me_yw_hw_cons);

		
		///////////////*				
		//Covariance between Consumption and Male Income
		
		local j=`j'+1;
		gen `eps_`j''	=	duc_duy - ( `k1'*`k5' * m_s_u1 * cov_dQ_bardQ  //Cov(q(1,t), qbar(1,t))
									+   `k2'*`k6' * m_s_u2 * cov_dQ_bardQ  //Cov(q(2,t), qbar(2,t))
									+   `k3'*`k7' * m_s_v1 * cov_dP_bardP  //Cov(p(1,t), pbar(1,t))
									+   `k4'*`k8' * m_s_v2 * cov_dP_bardP  //Cov(p(2,t), pbar(2,t))
									+ (`k1'*`k6'+`k2'*`k5') * m_r_u1u2 * cov_dQ_bardQ //Cov(q(1,t), qbar(2,t))
									+(`k3'*`k8'+`k4'*`k7') * m_r_v1v2 * cov_dP_bardP );   //Cov(p(1,t), pbar(2,t))  //no measurement error  Checked 10-7
						 
	// 		gen `eps_`j''	=	duc_duy - 2*(`k1'*`k5'*m_s_u1 + `k2'*`k6'*m_s_u2
	// 						+ (`k1'*`k6'+`k2'*`k5')*(m_r_u1u2)
	// 						+ `k3'*`k7'*m_s_v1 + `k4'*`k8'*m_s_v2
	// 						+ (`k3'*`k8' + `k4'*`k7')*(m_r_v1v2)); 
		
		///////////////*				
		//Covariance between Consumption and Female Income
		
		local j=`j'+1;
		gen `eps_`j''	=	duc_duyw - ( `k1'*`k9' * m_s_u1 * cov_dQ_bardQ  //Cov(q(1,t), qbar(1,t))
									+   `k2'*`k10' * m_s_u2 * cov_dQ_bardQ  //Cov((2,t), qbar(2,t))
									+   `k3'*`k11' * m_s_v1 * cov_dP_bardP  //Cov(p(1,t), pbar(1,t))
									+   `k4'*`k12' * m_s_v2 * cov_dP_bardP //Cov(p(2,t), pbar(2,t))
									+ (`k1'*`k10'+`k2'*`k9') * m_r_u1u2  * cov_dQ_bardQ  //Cov(q(1,t), qbar(2,t))
									+(`k3'*`k12'+`k4'*`k11') * m_r_v1v2  * cov_dP_bardP  );  //Cov(p(1,t), pbar(2,t))   Checked 10-7

	// 		gen `eps_`j''	=	duc_duyw - (`k1'*`k9'*(m_s_u1)*(2-(r_u1)^2*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
	// 						+ `k2'*`k10'*(m_s_u2)*(2-(r_u2)^2*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
	// 						+ (`k1'*`k10' + `k2'*`k9')*(2*m_r_u1u2 - (r_u1*r_u2)*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
	// 						+ `k3'*`k11'*m_s_v1*(2-(r_v1)^2*z_eta*inv_mills)
	// 						+ `k4'*`k12'*m_s_v2*(2-(r_v2)^2*z_eta*inv_mills)
	// 						+ (`k3'*`k12' + `k4'*`k11')*(2*m_r_v1v2 - (r_v1*r_v2)*z_eta*inv_mills));

		
		///////////////*				
		//Covariance between Male Income and Female Income
		local j=`j'+1;
		gen `eps_`j''	=	duy_duyw - ( `k5'*`k9' * m_s_u1 * var_bardQ  //Var(qbar(1,t)) 
									 +   `k6'*`k10' * m_s_u2 * var_bardQ //Var(qbar(2,t))
									 +   `k7'*`k11' * m_s_v1 * var_bardP //Var(pbar(1,t))
									 +   `k8'*`k12' * m_s_v2 * var_bardP //Var(pbar(2,t))
									 + (`k5'*`k10' + `k6'*`k9') * m_r_u1u2 * var_bardQ  //Cov(qbar(1,t),qbar(2,t))
									 + (`k7'*`k12' + `k8'*`k11') * m_r_v1v2 * var_bardP  ); //Cov(pbar(1,t), pbar(2,t))  //no measurement error  CHECKED 10-7

	// 		gen `eps_`j''	=	duy_duyw - (`k5'*`k9'*(m_s_u1)*(2-(r_u1)^2*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
	// 						+ `k6'*`k10'*(m_s_u2)*(2-(r_u2)^2*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
	// 						+ (`k5'*`k10' + `k6'*`k9')*(2*m_r_u1u2 - (r_u1*r_u2)*(z_eta*inv_mills + l2_z_eta*l2_inv_mills))
	// 						+ `k7'*`k11'*m_s_v1*(2-(r_v1)^2*z_eta*inv_mills)
	// 						+ `k8'*`k12'*m_s_v2*(2-(r_v2)^2*z_eta*inv_mills)
	// 						+ (`k7'*`k12' + `k8'*`k11')*(2*m_r_v1v2 - (r_v1*r_v2)*z_eta*inv_mills));

		///////////////*				END OF BATCH 1

		
		//Covariance between Consumption and lagged consumption
		local j=`j'+1;
		gen `eps_`j''	=	duc_lag - ( (`k1')^2 * l2_s_u1 * cov_dQ_LdQ  //Cov(q(1,t), q(1, t-2))
									+   (`k2')^2 * l2_s_u2 * cov_dQ_LdQ  //Cov(q(2,t), q(2,t-2))
									+ (2*`k1'*`k2') * l2_r_u1u2 * cov_dQ_LdQ //Cov(q(1,t), q(2,t-2))  //Checked 10-7
									- `s_omega'  );  //what is s_omega

 		//gen `eps_`j''	=	duc_lag - (-(`k1')^2*(l2_s_u1) -(`k2')^2*(l2_s_u2) -`s_omega' - 2*(`k1'*`k2')*(l2_r_u1u2));


		///////////////*
		//Covariance between Male Wage and lagged Consumption. 
		local j=`j'+1;
		gen `eps_`j''	=  	duw_duc_lag - ( `k1' * l2_s_u1 * cov_bardQ_dQL  //cov(qbar(1,t), q(1,t-2))
										+   `k2' * l2_r_u1u2 * cov_bardQ_dQL); // cov (qbar (1,t) , q(2, t-2))  CHECKED 10-7

	// 		gen `eps_`j''	=  	duw_duc_lag - (-`k1'*l2_s_u1
	// 						- `k2'*(l2_r_u1u2));


		///////////////*
		//Covariance between Female Wage and lagged Consumption
		local j=`j'+1;
		gen `eps_`j''	=  	duww_duc_lag - ( `k1' * l2_r_u1u2 * cov_bardQ_dQL    // cov(qbar(2,t), q(1, t-2))
										 +  `k2' * l2_s_u2 * cov_bardQ_dQL  ); // cov(qbar(2,t), q(2, t-2))   CHECKED 10-7
						
	// 		gen `eps_`j''	=  	duww_duc_lag - (-`k1'*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills)
	// 						-`k2'*(l2_s_u2)*(1-(r_u2)^2*l2_z_eta*l2_inv_mills));

	
		///////////////*	
		//Covariance between Consumption and lagged male wage
		local j=`j'+1;
		gen `eps_`j''	=  	duc_duw_lag - ( `k1' * l2_s_u1 * cov_bardQ_dQL      //cov(q(1,t), qbar(1,t-2))
										+	`k2' * l2_r_u1u2 * cov_bardQ_dQL ) ; //cov(q(2,t), qbar(1,t-2))  CHECKED 10-7

	// 		gen `eps_`j''	=  	duc_duw_lag - (-`k1'*l2_s_u1
	// 						- `k2'*(l2_r_u1u2));


		///////////////*
		//Covariance between Consumption and lagged female wage
		local j=`j'+1;
		gen `eps_`j''	=  	duc_duww_lag - ( `k1' * l2_r_u1u2 * cov_bardQ_dQL      //cov (q(1,t), qbar(2, t-2))
										 +	 `k2' * l2_s_u2 * cov_bardQ_dQL   ) ;  //cov (q(2,t), qbar(2, t-2))  CHECKED 10-7

	// 		local j=`j'+1;
	// 		gen `eps_`j''	=  	duc_duww_lag - (-`k1'*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills)
	// 						- `k2'*(l2_s_u2)*(1-(r_u2)^2*l2_z_eta*l2_inv_mills));

			//END OF BATCH 2	
		
		///////////////*		
		//Covariance between Male Income and Lagged Consumption
		local j=`j'+1;
		gen `eps_`j''	=	duy_duc_lag - ( (`k1'*`k5') * l2_s_u1 * cov_bardQ_dQL      // cov (qbar(1,t), q(1,t-2))
										+	(`k2'*`k6') * l2_s_u2 * cov_bardQ_dQL      // cov (qbar(2,t), q(2,t-2))
										+ (`k1'*`k6' + `k2'*`k5') * l2_r_u1u2 * cov_bardQ_dQL  ); // cov(qbar(1,t), q(2,t-2)) then again but flipped.   CHECKED 10-7
						
	// 		gen `eps_`j''	=	duy_duc_lag - (-`k1'*`k5'*l2_s_u1
	// 						- (`k1'*`k6'+`k2'*`k5')*(l2_r_u1u2)
	// 						- `k2'*`k6'*l2_s_u2);


		
		///////////////*
		//Covariance between Female Income and Lagged Consumption. 
		
		local j=`j'+1;
		gen `eps_`j''	=	duyw_duc_lag - ( (`k1'*`k9') * l2_s_u1 * cov_bardQ_dQL    // cov (qbar(1,t), q(1,t-2))
										 +   (`k2'*`k10') * l2_s_u2 * cov_bardQ_dQL 	// cov (qbar(2,t), q(2,t-2))
										 + (`k1'*`k10' + `k2'+`k9') * l2_r_u1u2 * cov_bardQ_dQL );  // cov(qbar(1,t), q(2,t-2)) then again but flipped.  CHECKED 10-7
						
	// 		gen `eps_`j''	=	duyw_duc_lag - (-`k1'*`k9'*l2_s_u1*(1-(r_u1)^2*l2_z_eta*l2_inv_mills)
	// 						- (`k1'*`k10'+`k2'*`k9')*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills)
	// 						- `k2'*`k10'*l2_s_u2*(1-(r_u2)^2*l2_z_eta*l2_inv_mills));

	
		///////////////*
		//Covariance between Consumption and Lagged Male Income 
		local j=`j'+1;

		gen `eps_`j''	=	duc_duy_lag - ( (`k1'*`k5') * l2_s_u1 * cov_bardQ_dQL     // cov (q(1,t), qbar(1,t-2))
										+   (`k2'*`k6') * l2_s_u2 * cov_bardQ_dQL     // cov (q(2,t), qbar(2,t-2))
										+ (`k1'*`k6' +`k2'*`k5') * l2_r_u1u2 * cov_bardQ_dQL );  // cov(q(1,t), qbar(2,t-2)) then again but flipped.  CHECKED 10-7
											
	// 		gen `eps_`j''	=	duc_duy_lag - (-`k1'*`k5'*l2_s_u1
	// 						- (`k1'*`k6'+`k2'*`k5')*(l2_r_u1u2)
	// 						- `k2'*`k6'*l2_s_u2);


		
		///////////////*
		//Covariance between Consumption and Lagged Female Income
		local j=`j'+1;
		gen `eps_`j''	=	duc_duyw_lag - ( (`k1'*`k9') * l2_s_u1 * cov_bardQ_dQL     	// cov (q(1,t), qbar(1,t-2))
										 +   (`k2'*`k10') * l2_s_u2 * cov_bardQ_dQL 	// cov (q(2,t), qbar(2,t-2))
										 + (`k1'*`k10' + `k2'*`k9') * l2_r_u1u2 * cov_bardQ_dQL   ); // cov(q(1,t), qbar(2,t-2)) then again but flipped.   CHECKED 10-7
						
	// 		gen `eps_`j''	=	duc_duyw_lag - (-`k1'*`k9'*l2_s_u1*(1-(r_u1)^2*l2_z_eta*l2_inv_mills)
	// 						- (`k1'*`k10'+`k2'*`k9')*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills)
	// 						- `k2'*`k10'*l2_s_u2*(1-(r_u2)^2*l2_z_eta*l2_inv_mills));

		
		///////////////*
		// Covariance between Wage and Lagged Male Income
		local j=`j'+1;
		gen `eps_`j''	=  	duw_duy_lag - ( `k5' * l2_s_u1 * cov_bardQ_LbardQ   // cov(qbar(1,t), qbar(1,t-2))
										+   `k6' * l2_r_u1u2 * cov_bardQ_LbardQ  ) ; // cov(qbar(1,t), qbar(2, t-2))  CHECKED 10-7
						
	// 		gen `eps_`j''	=  	duw_duy_lag - (-`k5'*(l2_s_u1)
	// 						- `k6'*(l2_r_u1u2)
	// 						- (s_me_y - cov_me_y_h));

					
		///////////////*
		//Covariance between Female Wage and Lagged Female Income
		local j=`j'+1;
		gen `eps_`j''	=  	duww_duyw_lag - ( `k9' * l2_r_u1u2 * cov_bardQ_LbardQ  // cov(qbar(2,t), qbar(1,t-2))
										  +	  `k10' * l2_s_u2 * cov_bardQ_LbardQ   ); // cov(qbar(2,t), qbar(2,t-2))  CHECKED 10-7

	// 		gen `eps_`j''	=  	duww_duyw_lag - (-`k9'*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills)
	// 						- `k10'*(l2_s_u2)*(1-(r_u2)^2*l2_z_eta*l2_inv_mills)
	// 						- (s_me_yw_cons - cov_me_yw_hw_cons));
						

		///////////////*
		//Covariance between Male Income and Lagged Male Wage
		local j=`j'+1;
		gen `eps_`j''	=  	duy_duw_lag - ( `k5' * l2_s_u1 * cov_bardQ_LbardQ  //cov(qbar(1,t), qbar(1,t-2))
										+   `k6' * l2_r_u1u2 * cov_bardQ_LbardQ  ); //cov(qbar(1,t), qbar(2,t-2))   CHECKED 10-7

	// 		gen `eps_`j''	=  	duy_duw_lag - (-`k5'*l2_s_u1
	// 						- `k6'*(l2_r_u1u2)
	// 						- (s_me_y - cov_me_y_h));

						
		///////////////*	
		//Covariance between Female Income and Lagged Female Wage
		local j=`j'+1;
		gen `eps_`j''	=  	duyw_duww_lag - ( `k9' * l2_r_u1u2 * cov_bardQ_LbardQ  //cov(qbar(1,t), qbar(2,t-2)) 
										  +   `k10' * l2_s_u2 * cov_bardQ_LbardQ  ); //cov(qbar(2,t), qbar(2,t-2))  CHECKED 10-7

	// 		gen `eps_`j''	=  	duyw_duww_lag - (-`k9'*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills)
	// 						-`k10'*(l2_s_u2)*(1-(r_u2)^2*l2_z_eta*l2_inv_mills)
	// 						- (s_me_yw_cons - cov_me_yw_hw_cons));
	
		//End of the Third Batch. 
						
		///////////////*
		//Covariance between Male Wage and Lagged Female Income
		local j=`j'+1;
		gen `eps_`j''	=  	duw_duyw_lag - ( `k9' * l2_s_u1 * cov_bardQ_LbardQ   //cov(qbar(1,t), qbar(1,t-2))
										 +   `k10' * l2_r_u1u2 * cov_bardQ_LbardQ   );  //cov(qbar(1,t), qbar(2,t-2))  CHECKED 10-7
						
	// 		gen `eps_`j''	=  	duw_duyw_lag - (-`k9'*(l2_s_u1)*(1-(r_u1)^2*l2_z_eta*l2_inv_mills)
	// 						-`k10'*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills));

						
		///////////////*				
		//Covariance between Female Wage and Lagged Male Income				
		local j=`j'+1;
		gen `eps_`j''	=  	duww_duy_lag - (`k5' * l2_r_u1u2 * cov_bardQ_LbardQ  //cov(qbar(2,t), qbar(1,t-2))
										 +  `k6' * l2_s_u2 * cov_bardQ_LbardQ  ); //cov(qbar(2,t), qbar(2,t-2))  CHECKED 10-7 

	// 		gen `eps_`j''	=  	duww_duy_lag - (-`k5'*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills)
	// 						-`k6'*(l2_s_u2)*(1-(r_u2)^2*l2_z_eta*l2_inv_mills));



		///////////////*
		//Covariance between Female Income and Lagged Male Wage
		local j=`j'+1;
		gen `eps_`j''	=  	duyw_duw_lag - ( `k9' * l2_s_u1 * cov_bardQ_LbardQ  //cov(qbar(1,t), qbar(1,t-2)) 
										 +   `k10' * l2_r_u1u2 * cov_bardQ_LbardQ  ); //cov(qbar(2,t), qbar(1,t-2))	CHECKED 10-7					
	// 		gen `eps_`j''	=  	duyw_duw_lag - (-`k9'*(l2_s_u1)*(1-(r_u1)^2*l2_z_eta*l2_inv_mills)
	// 						-`k10'*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills));

						
		///////////////*
		//Covariance between Male Income and Lagged Female Wage				
		local j=`j'+1;
		gen `eps_`j''	=  	duy_duww_lag - ( `k5' * l2_r_u1u2 * cov_bardQ_LbardQ  //cov(qbar(1,t), qbar(2,t-2))
										 +   `k6' * l2_s_u2 * cov_bardQ_LbardQ   ); //cov(qbar(2,t), qbar(2,t-2))     CHECKED 10-7
												
	// 		gen `eps_`j''	=  	duy_duww_lag - (-`k5'*(l2_r_u1u2 - r_u1*r_u2*l2_z_eta*l2_inv_mills)
	// 						-`k6'*(l2_s_u2)*(1-(r_u2)^2*l2_z_eta*l2_inv_mills));

						

						
						
		#delimit cr
		
		local i=1
		foreach var of varlist `varlist' {
			replace `var' = `eps_`i''
			local i=`i'+1
		}
	}
end

mata:
void calc_invQR()
{
	q= st_data((1,.),("l2_q"))
	tau= st_data((1,.),("tau"))

	vpsi_c_lam	= st_local("psi_c_lam")
	vpsi_c_w1	= st_local("psi_c_w1")
	vpsi_c_w2	= st_local("psi_c_w2")
	vpsi_h1_lam	= st_local("psi_h1_lam")
	vpsi_h1_w1	= st_local("psi_h1_w1")
	vpsi_h1_w2	= st_local("psi_h1_w2")
	vpsi_h2_lam	= st_local("psi_h2_lam")
	vpsi_h2_w1	= st_local("psi_h2_w1")
	vpsi_h2_w2	= st_local("psi_h2_w2")
	
	psi_c_lam	= st_data((1,.),(vpsi_c_lam))
	psi_c_w1	= st_data((1,.),(vpsi_c_w1))
	psi_c_w2	= st_data((1,.),(vpsi_c_w2))
	psi_h1_lam	= st_data((1,.),(vpsi_h1_lam))
	psi_h1_w1	= st_data((1,.),(vpsi_h1_w1))
	psi_h1_w2	= st_data((1,.),(vpsi_h1_w2))
	psi_h2_lam	= st_data((1,.),(vpsi_h2_lam))
	psi_h2_w1	= st_data((1,.),(vpsi_h2_w1))
	psi_h2_w2	= st_data((1,.),(vpsi_h2_w2))

	eta_c_p		= strtoreal(st_local("eta_c_p"))
	eta_c_w1	= strtoreal(st_local("eta_c_w1"))
	eta_h1_w1	= strtoreal(st_local("eta_h1_w1"))
	eta_c_w2	= strtoreal(st_local("eta_c_w2"))
	eta_h1_w2	= strtoreal(st_local("eta_h1_w2"))
	eta_h2_w2	= strtoreal(st_local("eta_h2_w2"))
	eta_h1_p	= strtoreal(st_local("eta_h1_p"))
	eta_h2_p	= strtoreal(st_local("eta_h2_p"))
	eta_h2_w1 	= strtoreal(st_local("eta_h2_w1"))

	Q=J(3,3,0)
	R=J(3,3,0)

	for (i=1 ; i<=rows(q) ; i++) {
		Q[1,1]= 1
		Q[1,2]= tau[i]*q[i]*(eta_c_w1 + eta_c_w2)
		Q[1,3]= tau[i]*(1-q[i])*(eta_c_w1 + eta_c_w2)
		Q[2,2]= 1+tau[i]*q[i]*(eta_h1_w1 + eta_h1_w2)
		Q[2,3]= tau[i]*(1-q[i])*(eta_h1_w1 + eta_h1_w2)
		Q[3,2]= tau[i]*q[i]*(eta_h2_w1 + eta_h2_w2)
		Q[3,3]= 1+tau[i]*(1-q[i])*(eta_h2_w1 + eta_h2_w2)

		R[1,1]= -eta_c_p + eta_c_w1 + eta_c_w2
		R[1,2]= eta_c_w1 - tau[i]*q[i]*(eta_c_w1 + eta_c_w2)
		R[1,3]= eta_c_w2 - tau[i]*(1-q[i])*(eta_c_w1 + eta_c_w2)
		R[2,1]= eta_h1_p + eta_h1_w1 + eta_h1_w2
		R[2,2]= eta_h1_w1 - tau[i]*q[i]*(eta_h1_w1 + eta_h1_w2)
		R[2,3]= eta_h1_w2 - tau[i]*(1-q[i])*(eta_h1_w1 + eta_h1_w2)
		R[3,1]= eta_h2_p + eta_h2_w1 + eta_h2_w2
		R[3,2]= eta_h2_w1 - tau[i]*q[i]*(eta_h2_w1 + eta_h2_w2)
		R[3,3]= eta_h2_w2 - tau[i]*(1-q[i])*(eta_h2_w1 + eta_h2_w2)
		
		PSI= pinv(Q)*R
		
		psi_c_lam[i]=PSI[1,1]
		psi_c_w1[i]=PSI[1,2]
		psi_c_w2[i]=PSI[1,3]
		psi_h1_lam[i]=PSI[2,1]
		psi_h1_w1[i]=PSI[2,2]
		psi_h1_w2[i]=PSI[2,3]
		psi_h2_lam[i]=PSI[3,1]
		psi_h2_w1[i]=PSI[3,2]
		psi_h2_w2[i]=PSI[3,3]

	}

	st_store(.,vpsi_c_lam,psi_c_lam)
	st_store(.,vpsi_c_w1,psi_c_w1)
	st_store(.,vpsi_c_w2,psi_c_w2)
	st_store(.,vpsi_h1_lam,psi_h1_lam)
	st_store(.,vpsi_h1_w1,psi_h1_w1)
	st_store(.,vpsi_h1_w2,psi_h1_w2)
	st_store(.,vpsi_h2_lam,psi_h2_lam)
	st_store(.,vpsi_h2_w1,psi_h2_w1)
	st_store(.,vpsi_h2_w2,psi_h2_w2)
}
end
