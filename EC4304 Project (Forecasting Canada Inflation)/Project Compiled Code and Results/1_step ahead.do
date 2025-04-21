//Part1
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"
tsmktim time, start(1981m1)
gen t=_n

//ac and pac plots
//ac ld1_cpi_all_can
//pac ld1_cpi_all_can

//1-step ahead
//AR model for 1 step ahead forecast
* Set the desired values for starting and ending t and AR range
local start_t = 15
local end_t = 268
local ar_start = 1
local ar_end = 14

* Create an empty local macro to later store the AR model names for the estimates stats command
local models ""

* Loop through the AR orders and run the regressions
forvalues i = `ar_start'/`ar_end' {
    qui reg ld1_cpi_all_can L(1/`i').ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
    estimates store AR`i'
    * Append the model name to the models macro
    local models "`models' AR`i'"
}

estimates stats `models'


//Best AR based on AIC: AR12, Best AR based on BIC: AR10


//Building ADL models: Granger Causality Test, Based on AIC for AR
disp 0.75*254^(1/3) /*SW lag default for HAC*/

forvalues i = 1/12 { 
qui newey ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/`i').m_base1 if t >= `start_t' & t <= `end_t', lag(4)
testparm L(1/`i').m_base1 
} 
//significant at 3, we can include up to 12 to test

forvalues i = 1/12 { 
qui newey ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/`i').m2p if t >= `start_t' & t <= `end_t', lag(4)
testparm L(1/`i').m2p 
} //significant up till 12, we test up to 12 lags of m2p

forvalues i = 1/12 { 
qui newey ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/`i').m3 if t >= `start_t' & t <= `end_t', lag(4)
testparm L(1/`i').m3
} //no significant lags of m3, we do not include

forvalues i = 1/12 { 
qui newey ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/`i').bank_rate_l if t >= `start_t' & t <= `end_t', lag(4)
testparm L(1/`i').bank_rate_l
} //bank rate l significant up to lag 5, we test till 12

forvalues i = 1/12 { 
qui newey ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/`i').tbill_6mbank_rate if t >= `start_t' & t <= `end_t', lag(4)
testparm L(1/`i').tbill_6mbank_rate
} //significant up till 12 we test up to 12 lags of tbill_6mbank_rate

forvalues i = 1/12 { 
qui newey ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/`i').ippi_can  if t >= `start_t' & t <= `end_t', lag(4)
testparm L(1/`i').ippi_can 
} //significnat up till 12, we test up to 12 lags of ippi_can

forvalues i = 1/12 { 
qui newey ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/`i').unemp_can  if t >= `start_t' & t <= `end_t', lag(4)
testparm L(1/`i').unemp_can 
} //significnat at lag 7, we test up to 12 lags of unemp_can


local models_mbase ""

//testing to find best m_base1
local ar_range = "12"
local adl_range = "1 2 3 4 5 6 7 8 9 10 11 12"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld1_cpi_all_can L(1/`ar').ld1_cpi_all_can L(1/`adl').m_base1 if t >= `start_t' & t <= `end_t', robust
        estimates store ADL_mbase_`adl'
        * Append the model name to the models macro
        local models_mbase "`models_mbase' ADL_mbase_`adl'"
    }
}

estimates stats AR12 `models_mbase' // mbase with 3 lags helps improve ar12 model slightly


local models_m2p ""

//testing to find best m2
local ar_range = "12"
local adl_range = "1 2 3 4 5 6 7 8 9 10 11 12"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld1_cpi_all_can L(1/`ar').ld1_cpi_all_can L(1/`adl').m2p if t >= `start_t' & t <= `end_t', robust
        estimates store ADL_m2p_`adl'
        * Append the model name to the models macro
        local models_m2p "`models_m2p' ADL_m2p_`adl'"
    }
}

estimates stats AR12 `models_m2p' //m2p 1 lags is the best


local models_br ""

//testing to find best br
local ar_range = "12"
local adl_range = "1 2 3 4 5 6 7 8 9 10 11 12"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld1_cpi_all_can L(1/`ar').ld1_cpi_all_can L(1/`adl').bank_rate_l if t >= `start_t' & t <= `end_t', robust
        estimates store ADL_br_`adl'
        * Append the model name to the models macro
        local models_br "`models_br' ADL_br_`adl'"
    }
}

estimates stats AR12 `models_br' // 1 lag of bank rate improves



local models_tbill ""

//testing to find best tbill
local ar_range = "12"
local adl_range = "1 2 3 4 5 6 7 8 9 10 11 12"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld1_cpi_all_can L(1/`ar').ld1_cpi_all_can L(1/`adl').tbill_6mbank_rate if t >= `start_t' & t <= `end_t', robust
        estimates store ADL_tbill_`adl'
        * Append the model name to the models macro
        local models_tbill "`models_tbill' ADL_tbill_`adl'"
    }
}

estimates stats AR12 `models_tbill' //tbill does not improve AR12


local models_unemp ""

//testing to find best unemp
local ar_range = "12"
local adl_range = "1 2 3 4 5 6 7 8 9 10 11 12"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld1_cpi_all_can L(1/`ar').ld1_cpi_all_can L(1/`adl').unemp_can if t >= `start_t' & t <= `end_t', robust
        estimates store ADL_unemp_`adl'
        * Append the model name to the models macro
        local models_unemp "`models_unemp' ADL_unemp_`adl'"
    }
}

estimates stats AR12 `models_unemp' // unemp lag 7 is the best


local models_ippi ""

//testing to find best m_base1
local ar_range = "12"
local adl_range = "1 2 3 4 5 6 7 8 9 10 11 12"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld1_cpi_all_can L(1/`ar').ld1_cpi_all_can L(1/`adl').ippi_can if t >= `start_t' & t <= `end_t', robust
        estimates store ADL_ippi_`adl'
        * Append the model name to the models macro
        local models_ippi "`models_ippi' ADL_ippi_`adl'"
    }
}

estimates stats AR12 `models_ippi' // 1 lag of ippi is the best

* Define the list of variables with their respective lags
local m1_chosen "L(1/3).m_base1"
local m2_chosen "L(1/2).m2p"
local br_chosen "L(1/1).bank_rate_l"
local un_chosen "L(1/7).unemp_can"
local ip_chosen "L(1/1).ippi_can"

* Test all combinations of 2 variables
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `m2_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1m2
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `br_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1br
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1un
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1ip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m2_chosen' `br_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2br
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m2_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2un
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m2_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2ip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brun
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `br_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `un_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store unip

* Test all combinations of 3 variables
* Test all combinations of 3 variables
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `m2_chosen' `br_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1m2br
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `m2_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1m2un
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `m2_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1m2ip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1brun
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `br_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1brip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `un_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1unip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m2_chosen' `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brun
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m2_chosen' `br_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m2_chosen' `un_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2unip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `br_chosen' `un_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brunip

* Test all combinations of 4 variables
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `m2_chosen' `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1m2brun
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `m2_chosen' `br_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1m2brip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `m2_chosen' `un_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1m2unip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `br_chosen' `un_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1brunip
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m2_chosen' `br_chosen' `un_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brunip

* Test the model with all 5 variables
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can `m1_chosen' `m2_chosen' `br_chosen' `un_chosen' `ip_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m1m2brunip

estimates stat AR1 AR12 ADL_mbase_3 ADL_m2p_1 ADL_br_1 ADL_unemp_7 ADL_ippi_1 ///
m1m2 m1br m1un m1ip m2br m2un m2ip brun brip ///
unip m1m2br m1m2un m1m2ip m1brun m1brip m1unip m2brun m2brip m2unip brunip ///
m1m2brun m1m2brip m1m2unip m1brunip m2brunip ///
m1m2brunip  //also check if same number of observations of all models


//best 3 ADL model is 1) m1m2brunip 2) m1m2unip 3) m1brunip 
estimates stats AR12 m1m2brunip m1m2unip m1brunip
//Part2
////////////////////////////////////////////////////////////////////////////////
//now we test using the validation set, get ready for forecast combination

//1-step ahead rolling window for AR models
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"

keep date ld1_cpi_all_can unemp_can ippi_can m2p bank_rate_l m_base1

tsmktim time, start(1981m1)
gen t=_n
local start_t = 15 //restrict start to make all models have same number of observations, for up to AR(14)
local end_t = 268
local validation_start = 269
local validation_end = 368
local windowsize = `validation_end' - `validation_start' //this is one less than the actual window size, to help with looping


//Generate random walk forecast
//note: for a given, rw is average of past h periods, so for 2 step ahead, forecast of t+3, is average of t, t-1, t-2
//please RW forecast accordingly
gen RW_h1 = L.ld1_cpi_all_can if t>= `validation_start' & t <=`validation_end' 


* Define the AR orders
local ar_orders = "1 12" //for both AR(1) and AR(12)
* Loop over each AR order
foreach ar_order in `ar_orders' {
    * Generate predictions for each AR order using a rolling window
	gen ar`ar_order'_h1 = .
    forvalues i = 0/`windowsize' {
        * Adjust the start and end for the training window
        local current_start_t = `start_t' + `i'
        local current_end_t = `end_t' + `i'
        
        * Estimate the AR model for the current training window
        qui reg ld1_cpi_all_can L(1/`ar_order').ld1_cpi_all_can if t >= `current_start_t' & t <= `current_end_t', robust
        
        * Predict the next value
        predict temp_pred, xb
        
        * Store the forecasted value in the pred_ar variable
        replace ar`ar_order'_h1 = temp_pred if t == `current_end_t' + 1
        
        * Drop the temporary prediction variable
        drop temp_pred
    }
}

//we also generate the forecasts for the models which are the best 3 adl models: 
gen m1m2brunip_h1= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i'
	
	* Estimate the AR model for the current training window
	qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/3).m_base1 L(1/2).m2p L(1/1).bank_rate_l L(1/8).unemp_can L(1/1).ippi_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m1m2brunip_h1 = temp_pred if t == `current_end_t' + 1
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 
gen m1m2unip_h1= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i'
	
	* Estimate the AR model for the current training window
	qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/3).m_base1 L(1/2).m2p L(1/8).unemp_can L(1/1).ippi_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m1m2unip_h1 = temp_pred if t == `current_end_t' + 1
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 
gen m1brunip_h1 = .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i'
	
	* Estimate the AR model for the current training window
	qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/3).m_base1 L(1/1).bank_rate_l L(1/8).unemp_can L(1/1).ippi_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m1brunip_h1 = temp_pred if t == `current_end_t' + 1
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//now we do forecast combination - Getting the weights
//We do forecast averaging, BMA and WAIC in the final step, because it doesnt need the real values.

//here we just do it for bates-granger and granger ramanathan as they need the true values
//Bates-granger
//first calculate the RMSE
gen true_values = ld1_cpi_all_can if t >= `validation_start' & t <= `validation_end'

egen e_RW_h1 = mean((RW_h1- true_values)^2)
gen rmse_RW_h1 = sqrt(e_RW_h1)

egen e_ar1_h1 = mean((ar1_h1- true_values)^2)
gen rmse_ar1_h1 = sqrt(e_ar1_h1)

egen e_ar12_h1 = mean((ar12_h1- true_values)^2)
gen rmse_ar12_h1 = sqrt(e_ar12_h1)

egen e_m1m2brunip_h1 = mean((m1m2brunip_h1- true_values)^2)
gen rmse_m1m2brunip_h1 = sqrt(e_m1m2brunip_h1)

egen e_m1m2unip_h1 = mean((m1m2unip_h1- true_values)^2)
gen rmse_m1m2unip_h1= sqrt(e_m1m2unip_h1)

egen e_m1brunip_h1 = mean((m1brunip_h1- true_values)^2)
gen rmse_m1brunip_h1 = sqrt(e_m1brunip_h1)

//Bates-Granger
//calculate variances of the different models to combine
//calculate weights of bates granger in the following way in excel
//1) for every forecast: calculate 1/MSE (not rmse)
//2) then calculate the weight by normalizing by sum across all forecasts
//code to output forecast mse and rmse into excel
scalar mse_RW = e_RW_h1[1]
scalar rmse_RW = rmse_RW_h1[1]

scalar mse_ar1 = e_ar1_h1[1]
scalar rmse_ar1 = rmse_ar1_h1[1]

scalar mse_ar12 = e_ar12_h1[1]
scalar rmse_ar12 = rmse_ar12_h1[1]

scalar mse_m1m2brunip= e_m1m2brunip_h1[1]
scalar rmse_m1m2brunip = rmse_m1m2brunip_h1[1]

scalar mse_m1m2unip = e_m1m2unip_h1[1]
scalar rmse_m1m2unip = rmse_m1m2unip_h1[1]

scalar mse_m1brunip = e_m1brunip_h1[1]
scalar rmse_m1brunip= rmse_m1brunip_h1[1]

* Create a numeric matrix with MSE and RMSE values
matrix F = (mse_RW, rmse_RW \ ///
		   mse_ar1, rmse_ar1 \ ///
           mse_ar12, rmse_ar12 \ ///
           mse_m1m2brunip, rmse_m1m2brunip \ /// 
           mse_m1m2unip, rmse_m1m2unip \ /// 
           mse_m1brunip, rmse_m1brunip) 
* Label the columns
matrix colnames F = mse rmse
* Label the rows with the forecast names
matrix rownames F = RW_h1 ar1_h1 ar12_h1 m1m2brunip m1m2unip m1brunip

* List the matrix to view it
matrix list F
* Export matrix F to Excel
putexcel set "Bates_Granger_weights_h1.xlsx", modify
putexcel A1=matrix(F), names

//weights for Bates-Granger
//RW = 0.105722425, AR1 = 0.178380597
//AR12 = 0.176931553, m1m2brunip = 0.180402181
//m1m2unip = 0.18277671, m1brunip = 0.176075615

//Granger-Ramanathan combination
*Unconstrained regression first to see which forecasts are collinear. If find any, drop them:
reg true_values RW_h1 ar1_h1 ar12_h1 m1m2brunip_h1 m1m2unip_h1 m1brunip_h1, noconstant
*Actually we find none, so we proceed to set the constraint, and use cnsreg to run the regression:
constraint 1 RW_h1 + ar1_h1 + ar12_h1 + m1m2brunip_h1 + m1m2unip_h1 + m1brunip_h1 = 1
cnsreg true_values RW_h1 ar1_h1 ar12_h1 m1m2brunip_h1 m1m2unip_h1 m1brunip_h1, constraints(1) noconstant
//dropping negative weights
constraint 2 RW_h1+ar1_h1+ar12_h1+m1m2unip_h1 = 1
cnsreg true_values RW_h1 ar1_h1 ar12_h1 m1m2unip_h1, constraints(2) noconstant
//dropping negative weights
constraint 3 RW_h1+ar1_h1+m1m2unip_h1 = 1
cnsreg true_values RW_h1 ar1_h1 m1m2unip_h1, constraints(3) noconstant

//Granger-Ramanathan weights
// RW_h1 = .027279, ar1_h1 = .4502477 , m1m2unip_h1 = .5224733

//We do forecast averaging, BMA and WAIC in the final step, because it doesnt need the real values.




//Part 3: same variable names, but now with all data up till before test period
////Final run of everything:
//1-step ahead rolling window for AR models
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"

keep date ld1_cpi_all_can unemp_can ippi_can m2p bank_rate_l m_base1

tsmktim time, start(1981m1)
gen t=_n
local start_t = 15 //restrict start to make all models have same number of observations, for up to AR(14)
local end_t = 368 //have all the data in for the test set
local test_start = 369
local test_end = 468
local windowsize = `test_end' - `test_start' //this is one less than the actual window size, to help with looping


//Generate random walk forecast'
//note: for a given, rw is average of past h periods, so for 2 step ahead, forecast of t+3, is average of t, t-1, t-2
//please RW forecast accordingly
gen RW_h1 = L.ld1_cpi_all_can if t>= `test_start' & t <=`test_end'  //note: for higher horize


* Define the AR orders
local ar_orders = "1 12" //for both AR(1) and AR(12)
* Loop over each AR order
foreach ar_order in `ar_orders' {
    * Generate predictions for each AR order using a rolling window
	gen ar`ar_order'_h1 = .
    forvalues i = 0/`windowsize' {
        * Adjust the start and end for the training window
        local current_start_t = `start_t' + `i'
        local current_end_t = `end_t' + `i'
        
        * Estimate the AR model for the current training window
        qui reg ld1_cpi_all_can L(1/`ar_order').ld1_cpi_all_can if t >= `current_start_t' & t <= `current_end_t', robust
        
        * Predict the next value
        predict temp_pred, xb
        
        * Store the forecasted value in the pred_ar variable
        replace ar`ar_order'_h1 = temp_pred if t == `current_end_t' + 1
        
        * Drop the temporary prediction variable
        drop temp_pred
    }
}

//we also generate the forecasts for the models which are the best 3 adl models: 
gen m1m2brunip_h1= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i'
	
	* Estimate the AR model for the current training window
	qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/3).m_base1 L(1/2).m2p L(1/1).bank_rate_l L(1/8).unemp_can L(1/1).ippi_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m1m2brunip_h1 = temp_pred if t == `current_end_t' + 1
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 
gen m1m2unip_h1= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i'
	
	* Estimate the AR model for the current training window
	qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/3).m_base1 L(1/2).m2p L(1/8).unemp_can L(1/1).ippi_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m1m2unip_h1 = temp_pred if t == `current_end_t' + 1
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models:
gen m1brunip_h1 = .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i'
	
	* Estimate the AR model for the current training window
	qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/3).m_base1 L(1/1).bank_rate_l L(1/8).unemp_can L(1/1).ippi_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m1brunip_h1 = temp_pred if t == `current_end_t' + 1
	
	* Drop the temporary prediction variable
	drop temp_pred
}

//generate forecast combinations:
//Forecast averaging:
gen forecast_average_h1 = (RW_h1 + ar1_h1 + ar12_h1 + m1m2brunip_h1 + m1m2unip_h1 + m1brunip_h1)/6

//Bates-Granger Combination:
//using weights from part 2
//RW = 0.105722425, AR1 = 0.178380597
//AR12 = 0.176931553, m1m2brunip = 0.180402181
//m1m2unip = 0.18277671, m1brunip = 0.176075615
gen bate_granger_h1 = 0.105722425*RW_h1 + 0.178380597*ar1_h1 + 0.176931553*ar12_h1 + 0.180402181*m1m2brunip_h1 + 0.18277671*m1m2unip_h1 + 0.176075615*m1brunip_h1

//Granger Ramanathan Combination:
//using weights from part 2
// RW_h1 = .027279, ar1_h1 = .4502477 , m1m2unip_h1 = .5224733
gen granger_raman_h1 = 0.027279*RW_h1 + 0.4502477*ar1_h1 + 0.5224733*m1m2unip_h1


//BMA of AR and ADL models
* Get the ICs of the different models
qui reg ld1_cpi_all_can L(1/1).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR1
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR12
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/3).m_base1 L(1/2).m2p L(1/1).bank_rate_l L(1/8).unemp_can L(1/1).ippi_can if t >= `start_t' & t <= `end_t', robust
estimates store m1m2brunip_h1
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/3).m_base1 L(1/2).m2p L(1/8).unemp_can L(1/1).ippi_can if t >= `start_t' & t <= `end_t', robust
estimates store m1m2unip_h1
qui reg ld1_cpi_all_can L(1/12).ld1_cpi_all_can L(1/3).m_base1 L(1/1).bank_rate_l L(1/8).unemp_can L(1/1).ippi_can if t >= `start_t' & t <= `end_t', robust
estimates store m1brunip_h1

estimates stats AR1 AR12 m1m2brunip_h1 m1m2unip_h1 m1brunip_h1

//weights following BMA. Copy paste BIC table in stata output into excel and calculate manually (select: right click copy table)
//weights according to BMA:
//AR1: 0.9987, AR12: 0.0013, others are too small, we leave output
gen bma_combined_h1 = 0.9987*ar1_h1 + 0.0013*ar12_h1


//WAIC of AR and ADL models
//weights following WAIC. Copy paste AIC table in stata output into excel and calculate manually (select: right click copy table)
//m1m2brunip = 0.3852177422, m1m2unip_h1 = 0.4462184984, m1brunip_h1 = 0.1685628604. others are too small
gen waic_combined_h1 = 0.3852177422* m1m2brunip_h1 + 0.4462184984*m1m2unip_h1 + 0.1685628604*m1brunip_h1


//now we do the DM test for each of the 10 forecasts against Random Walk
gen true_values = ld1_cpi_all_can if t >= `test_start' & t <= `test_end'

//errors of forecast
gen e_RW_h1 = true_values - RW_h1 
gen e_ar1_h1 = true_values - ar1_h1 
gen e_ar12_h1= true_values - ar12_h1 
gen e_m1m2brunip_h1 = true_values - m1m2brunip_h1 
gen e_m1m2unip_h1 = true_values - m1m2unip_h1 
gen e_m1brunip_h1 = true_values - m1brunip_h1 
gen e_forecast_average_h1 = true_values - forecast_average_h1 
gen e_bate_granger_h1 = true_values - bate_granger_h1 
gen e_granger_raman_h1 = true_values - granger_raman_h1 
gen e_bma_combined_h1= true_values - bma_combined_h1
gen e_waic_combined_h1= true_values - waic_combined_h1

//mse and rmse of errors
egen mse_RW_h1 = mean(e_RW_h1^2)
gen rmse_RW_h1 = sqrt(mse_RW_h1)

egen mse_ar1_h1 = mean(ar1_h1^2)
gen rmse_ar1_h1 = sqrt(mse_ar1_h1)

egen mse_ar12_h1 = mean(ar12_h1^2)
gen rmse_ar12_h1 = sqrt(mse_ar12_h1)

egen mse_m1m2brunip_h1 = mean(m1m2brunip_h1^2)
gen rmse_m1m2brunip_h1 = sqrt(mse_m1m2brunip_h1)

egen mse_m1m2unip_h1 = mean(m1m2unip_h1^2)
gen rmse_m1m2unip_h1 = sqrt(mse_m1m2unip_h1)

egen mse_m1brunip_h1 = mean(m1brunip_h1^2)
gen rmse_m1brunip_h1 = sqrt(mse_m1brunip_h1)

egen mse_forecast_average_h1 = mean(forecast_average_h1^2)
gen rmse_forecast_average_h1 = sqrt(mse_forecast_average_h1)

egen mse_bate_granger_h1 = mean(bate_granger_h1^2)
gen rmse_bate_granger_h1 = sqrt(mse_bate_granger_h1)

egen mse_granger_raman_h1 = mean(granger_raman_h1^2)
gen rmse_granger_raman_h1 = sqrt(mse_granger_raman_h1)

egen mse_bma_combined_h1 = mean(bma_combined_h1^2)
gen rmse_bma_combined_h1 = sqrt(mse_bma_combined_h1)

egen mse_waic_combined_h1 = mean(waic_combined_h1^2)
gen rmse_waic_combined_h1 = sqrt(mse_waic_combined_h1)


//getting loss differentials and the maxlag to use in DM test
gen d_rw_ar1 =e_RW_h1^2 - e_ar1_h1^2
ac d_rw_ar1
corrgram  d_rw_ar1
//RW-AR1 should use lag 0, looks like white noise
gen d_rw_ar12 =e_RW_h1^2 - e_ar12_h1^2
ac d_rw_ar12
corrgram  d_rw_ar12
//RW-AR12 should use lag 0, looks like white noise
gen d_rw_adl_1 =e_RW_h1^2 - e_m1m2brunip_h1^2
ac d_rw_adl_1
corrgram  d_rw_adl_1
//RW-m1m2brunip_h1 should use lag 0, looks like white noise
gen d_rw_adl_2 =e_RW_h1^2 - e_m1m2unip_h1^2
ac d_rw_adl_2
corrgram  d_rw_adl_2
//RW-m1m2unip_h1 should use lag 0, looks like white noise
gen d_rw_adl_3 =e_RW_h1^2 - e_m1brunip_h1^2
ac d_rw_adl_3
corrgram  d_rw_adl_3
//RW-m1brunip_h1 should use lag 0, looks like white noise
gen d_rw_fa =e_RW_h1^2 - e_forecast_average_h1^2
ac d_rw_fa
corrgram  d_rw_fa
//RW-forecast averaging should use lag 0, looks like white noise
gen d_rw_bg =e_RW_h1^2 - e_bate_granger_h1^2
ac d_rw_bg
corrgram  d_rw_bg
//RW-bates granger should use lag 0, looks like white noise
gen d_rw_gr =e_RW_h1^2 - e_granger_raman_h1^2
ac d_rw_gr
corrgram  d_rw_gr
//RW-granger ramanathan should use lag 0, looks like white noise
gen d_rw_bma =e_RW_h1^2 - e_bma_combined_h1^2
ac d_rw_bma
corrgram  d_rw_bma
//RW-bates granger should use lag 0, looks like white noise
gen d_rw_waic =e_RW_h1^2 - e_waic_combined_h1^2
ac d_rw_waic
corrgram  d_rw_waic
//RW-granger ramanathan should use lag 0, looks like white noise


//dm test against AR1
gen d_ar1_ar12 =e_ar1_h1^2 - e_ar12_h1^2
ac d_ar1_ar12
corrgram  d_ar1_ar12
//AR1-AR12 should use lag 0, looks like white noise
gen d_ar1_adl_1 =e_ar1_h1^2 - e_m1m2brunip_h1^2
ac d_ar1_adl_1
corrgram  d_ar1_adl_1
//AR1-m1m2brunip_h1 should use lag 0, looks like white noise
gen d_ar1_adl_2 =e_ar1_h1^2 - e_m1m2unip_h1^2
ac d_ar1_adl_2
corrgram  d_ar1_adl_2
//AR1-m1m2unip_h1 should use lag 0, looks like white noise
gen d_ar1_adl_3 =e_ar1_h1^2 - e_m1brunip_h1^2
ac d_ar1_adl_3
corrgram  d_ar1_adl_3
//AR1-m1brunip_h1 should use lag 0, looks like white noise
gen d_ar1_fa =e_ar1_h1^2 - e_forecast_average_h1^2
ac d_ar1_fa
corrgram  d_ar1_fa
//AR1-forecast averaging should use lag 0, looks like white noise
gen d_ar1_bg =e_ar1_h1^2 - e_bate_granger_h1^2
ac d_ar1_bg
corrgram  d_ar1_bg
//AR1-bates granger should use lag 0, looks like white noise
gen d_ar1_gr =e_ar1_h1^2 - e_granger_raman_h1^2
ac d_ar1_gr
corrgram  d_ar1_gr
//AR1-granger ramanathan should use lag 0, looks like white noise
gen d_ar1_bma =e_ar1_h1^2 - e_bma_combined_h1^2
ac d_ar1_bma
corrgram  d_ar1_bma
//AR1-bates granger should use lag 0, looks like white noise
gen d_ar1_waic =e_ar1_h1^2 - e_waic_combined_h1^2
ac d_ar1_waic
corrgram  d_ar1_waic
//AR1-granger ramanathan should use lag 0, looks like white noise



//compiling results with dm test:
matrix rmse_matrix = J(11 , 3, 0) //10 rows with 3 columns

//compiling rmse (make sure its in order)
matrix rmse_matrix[1, 1] = rmse_RW_h1[1]
matrix rmse_matrix[2, 1] = rmse_ar1_h1[1]
matrix rmse_matrix[3, 1] = rmse_ar12_h1[1]
matrix rmse_matrix[4, 1] = rmse_m1m2brunip_h1[1]
matrix rmse_matrix[5, 1] = rmse_m1m2unip_h1[1]
matrix rmse_matrix[6, 1] = rmse_m1brunip_h1[1]
matrix rmse_matrix[7, 1] = rmse_forecast_average_h1[1]
matrix rmse_matrix[8, 1] = rmse_bate_granger_h1[1]
matrix rmse_matrix[9, 1] = rmse_granger_raman_h1[1]
matrix rmse_matrix[10, 1] = rmse_bma_combined_h1[1]
matrix rmse_matrix[11, 1] = rmse_waic_combined_h1[1]

//DM test against RW (make sure its in order):
qui dmariano true_values RW_h1 ar1_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[2, 2] = r(p)
qui dmariano true_values RW_h1 ar12_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[3, 2] = r(p)
qui dmariano true_values RW_h1 m1m2brunip_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[4, 2] = r(p)
qui dmariano true_values RW_h1 m1m2unip_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[5, 2] = r(p)
qui dmariano true_values RW_h1 m1brunip_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[6, 2] = r(p)
qui dmariano true_values RW_h1 forecast_average_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[7, 2] = r(p)
qui dmariano true_values RW_h1 bate_granger_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[8, 2] = r(p)
qui dmariano true_values RW_h1 granger_raman_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[9, 2] = r(p)
qui dmariano true_values RW_h1 bma_combined_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[10, 2] = r(p)
qui dmariano true_values RW_h1 waic_combined_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[11, 2] = r(p)

//DM test against AR(1) (make sure its in order):
qui dmariano true_values ar1_h1 ar12_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[3, 3] = r(p)
qui dmariano true_values ar1_h1 m1m2brunip_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[4, 3] = r(p)
qui dmariano true_values ar1_h1 m1m2unip_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[5, 3] = r(p)
qui dmariano true_values ar1_h1 m1brunip_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[6, 3] = r(p)
qui dmariano true_values ar1_h1 forecast_average_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[7, 3] = r(p)
qui dmariano true_values ar1_h1 bate_granger_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[8, 3] = r(p)
qui dmariano true_values ar1_h1 granger_raman_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[9, 3] = r(p)
qui dmariano true_values ar1_h1 bma_combined_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[10, 3] = r(p)
qui dmariano true_values ar1_h1 waic_combined_h1, crit(mse) kernel(bartlett) maxlag(0)
matrix rmse_matrix[11, 3] = r(p)

* Name the columns and rows for clarity 
matrix colnames rmse_matrix = "RMSE" "DM RW" "DM AR(1)"
matrix rownames rmse_matrix =  RW_h1 AR1_h1 AR12_h1 ADL1_h1 ADL2_h1 ADL3_h1 FA_h1 BG_h1 GR_h1 BMA_h1 WAIC_h1

matrix list rmse_matrix

putexcel set "1_step_results.xlsx", replace
putexcel A1=matrix(rmse_matrix), names










