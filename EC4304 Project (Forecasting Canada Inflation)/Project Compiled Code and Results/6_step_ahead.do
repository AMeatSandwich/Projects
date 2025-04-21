//Part1
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"
tsmktim time, start(1981m1)
gen t=_n

//ac and pac plots
//ac ld6_cpi_all_can
//pac ld6_cpi_all_can

//1-step ahead
//AR model for 1 step ahead forecast
* Set the desired values for starting and ending t and AR range
local start_t = 25
local end_t = 268
local ar_start = 1
local ar_end = 17 //from ac graph

* Create an empty local macro to later store the AR model names for the estimates stats command
local models ""

* Loop through the AR orders and run the regressions
forvalues i = `ar_start'/`ar_end' {
    local x = `i' + 5
    qui reg ld6_cpi_all_can L(6/`x').ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
    estimates store AR`i'
    * Append the model name to the models macro
    local models "`models' AR`i'"
}

estimates stats `models'


//Best AR based on AIC: AR12, Best AR based on BIC: AR1


//Building ADL models: Granger Causality Test, Based on AIC for AR
disp 0.75*244^(1/3) /*SW lag default for HAC*/

local start_t = 25
local end_t = 268

forvalues i = 1/18 {
    local x = `i' + 5
	qui newey ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/`x').m_base1 if t >= `start_t' & t <= `end_t', lag(4)
	testparm L(6/`x').m_base1 
} 
//no significant lags of m_base1, we do not include

forvalues i = 1/18 {
    local x = `i' + 5 
	qui newey ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/`x').m2p if t >= `start_t' & t <= `end_t', lag(4)
	testparm L(6/`x').m2p 
} //significant up till 18, we test up to 18 lags of m2p


forvalues i = 1/18 { 
    local x = `i' + 5
	qui newey ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/`x').m3 if t >= `start_t' & t <= `end_t', lag(4)
	testparm L(6/`x').m3
} //significant up till lag 9, we test up to 12 lags of m3

forvalues i = 1/18 { 
    local x = `i' + 5
	qui newey ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/`x').bank_rate_l if t >= `start_t' & t <= `end_t', lag(4)
	testparm L(6/`x').bank_rate_l
} //no significant lags of bank rate l, we do not include

forvalues i = 1/18 { 
    local x = `i' + 5
	qui newey ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/`x').tbill_6mbank_rate if t >= `start_t' & t <= `end_t', lag(4)
	testparm L(6/`x').tbill_6mbank_rate
} //no significant lags of tbill_6mbank_rate, we do not include 

forvalues i = 1/18 { 
    local x = `i' + 5
	qui newey ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/`x').ippi_can  if t >= `start_t' & t <= `end_t', lag(4)
	testparm L(6/`x').ippi_can 
} //no significant lags of ippi_can, we do not include 

forvalues i = 1/18 { 
    local x = `i' + 5
	qui newey ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/`x').unemp_can  if t >= `start_t' & t <= `end_t', lag(4)
	testparm L(6/`x').unemp_can 
} //no significant lags of unemp_can, we do not include 

local start_t = 25
local end_t = 268

local models_m2p ""

//testing to find best m2p
local ar_range = "12"
local adl_range = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
    * Loop through the orders
    foreach adl in `adl_range' {
	    local x = `ar' + 5
		local y = `adl' + 5
        qui reg ld6_cpi_all_can L(6/`x').ld1_cpi_all_can L(6/`y').m2p if t >= `start_t' & t <= `end_t', robust
        estimates store ADL_m2p_`adl'
        * Append the model name to the models macro
        local models_m2p "`models_m2p' ADL_m2p_`adl'"
    }
}

estimates stats AR1 AR12 `models_m2p' //m2p 9 lags is the best


local models_br ""

//testing to find best m3
local ar_range = "12"
local adl_range = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
    * Loop through the orders
    foreach adl in `adl_range' {
	    local x = `ar' + 5
		local y = `adl' + 5
        qui reg ld6_cpi_all_can L(6/`x').ld1_cpi_all_can L(6/`y').m3 if t >= `start_t' & t <= `end_t', robust
        estimates store ADL_m3_`adl'
        * Append the model name to the models macro
        local models_m3 "`models_m3' ADL_m3_`adl'"
    }
}

estimates stats AR1 AR12 `models_m3' // 6 lag of m3 improves 

local start_t = 25
local end_t = 268



* Define the list of variables with their respective lags
local m3_chosen "L(6/11).m3"
local m2_chosen "L(6/14).m2p"

local start_t = 25
local end_t = 268

* Test all combinations of 2 variables
qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can `m3_chosen' `m2_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m3m2


estimates stat AR1 AR12 ADL_m3_6 ADL_m2p_9 m3m2  //also check if same number of observations of all models

//best 3 ADL model is 1) ADL_m2p_9 2) m3m2 3) ADL_m3_6 
estimates stat AR12 ADL_m2p_9 m3m2 ADL_m3_6
//Part2
////////////////////////////////////////////////////////////////////////////////
//now we test using the validation set, get ready for forecast combination

//1-step ahead rolling window for AR models
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"

keep date ld1_cpi_all_can ld6_cpi_all_can m2p m3

tsmktim time, start(1981m1)
gen t=_n

local start_t = 25 //restrict start to make all models have same number of observations, for up to AR(12)
local end_t = 268
local validation_start = 269
local validation_end = 368
local windowsize = `validation_end' - `validation_start' //this is one less than the actual window size, to help with looping


//Generate random walk forecast
//note: for a given, rw is average of past h periods, so for 2 step ahead, forecast of t+3, is average of t, t-1, t-2
//please RW forecast accordingly
cap gen RW_h6 = (L6.ld6_cpi_all_can + L7.ld6_cpi_all_can + L8.ld6_cpi_all_can + L9.ld6_cpi_all_can + L10.ld6_cpi_all_can + L11.ld6_cpi_all_can)/6 if t>= `validation_start' & t <=`validation_end' 


* Define the AR orders
local ar_orders = "1 12" //for both AR(1) and AR(12)
* Loop over each AR order
foreach ar_order in `ar_orders' {
    * Generate predictions for each AR order using a rolling window
	gen ar`ar_order'_h6 = .
    forvalues i = 0/`windowsize' {
        * Adjust the start and end for the training window
        local current_start_t = `start_t' + `i'
        local current_end_t = `end_t' + `i'- 6
        
        * Estimate the AR model for the current training window
		local x = `ar_order' + 5
        qui reg ld6_cpi_all_can L(6/`x').ld1_cpi_all_can if t >= `current_start_t' & t <= `current_end_t', robust
        
        * Predict the next value
        predict temp_pred, xb
        
        * Store the forecasted value in the pred_ar variable
        replace ar`ar_order'_h6 = temp_pred if t == `current_end_t' + 1 + 6
        
        * Drop the temporary prediction variable
        drop temp_pred
    }
}

//we also generate the forecasts for the models which are the best 3 adl models: 
gen m2_h6= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' -6
	
	* Estimate the AR model for the current training window
	qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/14).m2p if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2_h6 = temp_pred if t == `current_end_t' + 1 +6
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 
gen m3m2_h6= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' -6
	
	* Estimate the AR model for the current training window
	qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/14).m2p L(6/11).m3 if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m3m2_h6 = temp_pred if t == `current_end_t' + 1 + 6
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 
gen m3_h6 = .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' - 6
	
	* Estimate the AR model for the current training window
	qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/11).m3 if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m3_h6 = temp_pred if t == `current_end_t' + 1 + 6
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//now we do forecast combination - Getting the weights
//We do forecast averaging, BMA and WAIC in the final step, because it doesnt need the real values.

//here we just do it for bates-granger and granger ramanathan as they need the true values
//Bates-granger
//first calculate the RMSE
gen true_values = ld6_cpi_all_can if t >= `validation_start' & t <= `validation_end'

egen e_RW_h6 = mean((RW_h6- true_values)^2)
gen rmse_RW_h6 = sqrt(e_RW_h6)

egen e_ar1_h6 = mean((ar1_h6- true_values)^2)
gen rmse_ar1_h6 = sqrt(e_ar1_h6)

egen e_ar12_h6 = mean((ar12_h6- true_values)^2)
gen rmse_ar12_h6 = sqrt(e_ar12_h6)

egen e_m2_h6 = mean((m2_h6- true_values)^2)
gen rmse_m2_h6 = sqrt(e_m2_h6)

egen e_m3m2_h6 = mean((m3m2_h6- true_values)^2)
gen rmse_m3m2_h6= sqrt(e_m3m2_h6)

egen e_m3_h6 = mean((m3_h6- true_values)^2)
gen rmse_m3_h6 = sqrt(e_m3_h6)

//Bates-Granger
//calculate variances of the different models to combine
//calculate weights of bates granger in the following way in excel
//1) for every forecast: calculate 1/MSE (not rmse)
//2) then calculate the weight by normalizing by sum across all forecasts
//code to output forecast mse and rmse into excel
scalar mse_RW = e_RW_h6[1]
scalar rmse_RW = rmse_RW_h6[1]

scalar mse_ar1 = e_ar1_h6[1]
scalar rmse_ar1 = rmse_ar1_h6[1]

scalar mse_ar12 = e_ar12_h6[1]
scalar rmse_ar12 = rmse_ar12_h6[1]

scalar mse_m2= e_m2_h6[1]
scalar rmse_m2 = rmse_m2_h6[1]

scalar mse_m3m2 = e_m3m2_h6[1]
scalar rmse_m3m2 = rmse_m3m2_h6[1]

scalar mse_m3 = e_m3_h6[1]
scalar rmse_m3 = rmse_m3_h6[1]

* Create a numeric matrix with MSE and RMSE values
matrix F = (mse_RW, rmse_RW \ ///
		   mse_ar1, rmse_ar1 \ ///
           mse_ar12, rmse_ar12 \ ///
           mse_m2, rmse_m2 \ /// 
           mse_m3m2, rmse_m3m2 \ /// 
           mse_m3, rmse_m3) 
* Label the columns
matrix colnames F = mse rmse
* Label the rows with the forecast names
matrix rownames F = RW_h6 ar1_h6 ar12_h6 m2 m3m2 m3

* List the matrix to view it
matrix list F
* Export matrix F to Excel
putexcel set "Bates_Granger_weights_h6.xlsx", modify
putexcel A1=matrix(F), names

//weights for Bates-Granger
//RW = 0.118057256, AR1 = 0.198417061
//AR12 = 0.212672146, m2 = 0.146161821
//m3m2 = 0.149212056, m3 = 0.17547966

//Granger-Ramanathan combination
*Unconstrained regression first to see which forecasts are collinear. If find any, drop them:
reg true_values RW_h6 ar1_h6 ar12_h6 m2_h6 m3m2_h6 m3_h6, noconstant
*Actually we find none, so we proceed to set the constraint, and use cnsreg to run the regression:
constraint 1 RW_h6 + ar1_h6 + ar12_h6 + m2_h6 + m3m2_h6 + m3_h6 = 1
cnsreg true_values RW_h6 ar1_h6 ar12_h6 m2_h6 m3m2_h6 m3_h6, constraints(1) noconstant
//dropping negative weights
constraint 2 ar1_h6 + ar12_h6 + m2_h6 + m3_h6 = 1
cnsreg true_values ar1_h6 ar12_h6 m2_h6 m3_h6, constraints(2) noconstant
//dropping negative weights
constraint 3 ar1_h6 + ar12_h6 + m3_h6 = 1
cnsreg true_values ar1_h6 ar12_h6 m3_h6, constraints(3) noconstant

//Granger-Ramanathan weights
// ar1_h6 = .1968342  , ar12_h6 = .7152803 , m3_h6 = .0878855 

//We do forecast averaging, BMA and WAIC in the final step, because it doesnt need the real values.




//Part 3: same variable names, but now with all data up till before test period
////Final run of everything:
//1-step ahead rolling window for AR models
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"

keep date ld6_cpi_all_can ld1_cpi_all_can m2p m3

tsmktim time, start(1981m1)
gen t=_n

local start_t = 25 //restrict start to make all models have same number of observations, for up to AR(14)
local end_t = 368 //have all the data in for the test set
local test_start = 369
local test_end = 468
local windowsize = `test_end' - `test_start' //this is one less than the actual window size, to help with looping


//Generate random walk forecast'
//note: for a given, rw is average of past h periods, so for 2 step ahead, forecast of t+3, is average of t, t-1, t-2
//please RW forecast accordingly
gen RW_h6 = (L6.ld6_cpi_all_can + L7.ld6_cpi_all_can + L8.ld6_cpi_all_can + L9.ld6_cpi_all_can + L10.ld6_cpi_all_can + L11.ld6_cpi_all_can)/6 if t>= `test_start' & t <=`test_end'   //note: for higher horize


* Define the AR orders
local ar_orders = "1 12" //for both AR(1) and AR(12)
* Loop over each AR order
foreach ar_order in `ar_orders' {
    * Generate predictions for each AR order using a rolling window
	gen ar`ar_order'_h6 = .
    forvalues i = 0/`windowsize' {
        * Adjust the start and end for the training window
        local current_start_t = `start_t' + `i'
        local current_end_t = `end_t' + `i' - 6
        
        * Estimate the AR model for the current training window
		local x = `ar_order' + 5
        qui reg ld6_cpi_all_can L(6/`x').ld1_cpi_all_can if t >= `current_start_t' & t <= `current_end_t', robust
        
        * Predict the next value
        predict temp_pred, xb
        
        * Store the forecasted value in the pred_ar variable
        replace ar`ar_order'_h6 = temp_pred if t == `current_end_t' + 1 + 6
        
        * Drop the temporary prediction variable
        drop temp_pred
    }
}

//we also generate the forecasts for the models which are the best 3 adl models: 
gen m2_h6= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' - 6
	
	* Estimate the AR model for the current training window
	qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/14).m2p if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2_h6 = temp_pred if t == `current_end_t' + 1 +6
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 
gen m3m2_h6= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' - 6
	
	* Estimate the AR model for the current training window
	qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/14).m2p L(6/11).m3 if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m3m2_h6 = temp_pred if t == `current_end_t' + 1 + 6
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 
gen m3_h6 = .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' - 6
	
	* Estimate the AR model for the current training window
	qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/11).m3 if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m3_h6 = temp_pred if t == `current_end_t' + 1 + 6
	
	* Drop the temporary prediction variable
	drop temp_pred
}

//generate forecast combinations:
//Forecast averaging:
gen forecast_average_h6 = (RW_h6 + ar1_h6 + ar12_h6 + m2_h6 + m3m2_h6 + m3_h6)/6

//Bates-Granger Combination:
//using weights from part 2
//weights for Bates-Granger
//RW = 0.118057256, AR1 = 0.198417061
//AR12 = 0.212672146, m2 = 0.146161821
//m3m2 = 0.149212056, m3 = 0.17547966

gen bate_granger_h6 = 0.118057256*RW_h6 + 0.198417061*ar1_h6 + 0.212672146*ar12_h6 + 0.146161821*m2_h6 + 0.149212056*m3m2_h6 + 0.17547966*m3_h6

//Granger Ramanathan Combination:
//using weights from part 2
// ar1_h6 = .1968342  , ar12_h6 = .7152803 , m3_h6 = .0878855 

gen granger_raman_h6 = 0.1968342*ar1_h6 + 0.7152803*ar12_h6 + 0.0878855*m3_h6


//BMA of AR and ADL models
local start_t = 25 
local end_t = 368 -6 

* Get the ICs of the different models
qui reg ld6_cpi_all_can L(6/6).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR1
qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR12
qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/17).m2p if t >= `start_t' & t <= `end_t', robust
estimates store m2_h6
qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/17).m2p L(6/11).m3 if t >= `start_t' & t <= `end_t', robust
estimates store m3m2_h6
qui reg ld6_cpi_all_can L(6/17).ld1_cpi_all_can L(6/11).m3 if t >= `start_t' & t <= `end_t', robust
estimates store m3_h6

estimates stats AR1 AR12 m2_h6 m3m2_h6 m3_h6

//weights following BMA. Copy paste BIC table in stata output into excel and calculate manually (select: right click copy table)
//weights according to BMA:
//ar1_h6: 0.0052649345, ar12_h6: 0.9800158188 m3_h6: 0.0147179629, others are too small, we leave output
gen bma_combined_h6 = 0.002686913*ar1_h6 + 0.979353951*ar12_h6 + 0.017937493*m3_h6


//WAIC of AR and ADL models
//weights following WAIC. Copy paste AIC table in stata output into excel and calculate manually (select: right click copy table)
//m2_h6: 0.7060216940  m3m2_h6: 0.2137155585  m3_h6: 0.0802098211 others are too small
gen waic_combined_h6 = 0.838933687* m2_h6 + 0.15379671*m3m2_h6 + 0.00726546*m3_h6


//now we do the DM test for each of the 10 forecasts against Random Walk
gen true_values = ld6_cpi_all_can if t >= `test_start' & t <= `test_end'

//errors of forecast
gen e_RW_h6 = true_values - RW_h6 
gen e_ar1_h6 = true_values - ar1_h6
gen e_ar12_h6= true_values - ar12_h6
gen e_m2_h6 = true_values - m2_h6 
gen e_m3m2_h6 = true_values - m3m2_h6 
gen e_m3_h6 = true_values - m3_h6 
gen e_forecast_average_h6 = true_values - forecast_average_h6 
gen e_bate_granger_h6 = true_values - bate_granger_h6
gen e_granger_raman_h6 = true_values - granger_raman_h6 
gen e_bma_combined_h6= true_values - bma_combined_h6
gen e_waic_combined_h6= true_values - waic_combined_h6

//mse and rmse of errors
egen mse_RW_h6 = mean(e_RW_h6^2)
gen rmse_RW_h6 = sqrt(mse_RW_h6)

egen mse_ar1_h6 = mean(e_ar1_h6^2)
gen rmse_ar1_h6 = sqrt(mse_ar1_h6)

egen mse_ar12_h6 = mean(e_ar12_h6^2)
gen rmse_ar12_h6 = sqrt(mse_ar12_h6)

egen mse_m2_h6 = mean(e_m2_h6^2)
gen rmse_m2_h6 = sqrt(mse_m2_h6)

egen mse_m3m2_h6 = mean(e_m3m2_h6^2)
gen rmse_m3m2_h6 = sqrt(mse_m3m2_h6)

egen mse_m3_h6 = mean(e_m3_h6^2)
gen rmse_m3_h6 = sqrt(mse_m3_h6)

egen mse_forecast_average_h6 = mean(e_forecast_average_h6^2)
gen rmse_forecast_average_h6 = sqrt(mse_forecast_average_h6)

egen mse_bate_granger_h6 = mean(e_bate_granger_h6^2)
gen rmse_bate_granger_h6 = sqrt(mse_bate_granger_h6)

egen mse_granger_raman_h6 = mean(e_granger_raman_h6^2)
gen rmse_granger_raman_h6 = sqrt(mse_granger_raman_h6)

egen mse_bma_combined_h6 = mean(e_bma_combined_h6^2)
gen rmse_bma_combined_h6 = sqrt(mse_bma_combined_h6)

egen mse_waic_combined_h6 = mean(e_waic_combined_h6^2)
gen rmse_waic_combined_h6 = sqrt(mse_waic_combined_h6)


//getting loss differentials and the maxlag to use in DM test
gen d_rw_ar1 =e_RW_h6^2 - e_ar1_h6^2
ac d_rw_ar1
corrgram  d_rw_ar1
//RW-AR1 should use lag 1
gen d_rw_ar12 =e_RW_h6^2 - e_ar12_h6^2
ac d_rw_ar12
corrgram  d_rw_ar12
//RW-AR12 should use lag 7
gen d_rw_adl_1 =e_RW_h6^2 - e_m2_h6^2
ac d_rw_adl_1
corrgram  d_rw_adl_1
//RW-m2_h6 should use lag 7
gen d_rw_adl_2 =e_RW_h6^2 - e_m3m2_h6^2
ac d_rw_adl_2
corrgram  d_rw_adl_2
//RW-m3m2_h6 should use lag 7
gen d_rw_adl_3 =e_RW_h6^2 - e_m3_h6^2
ac d_rw_adl_3
corrgram  d_rw_adl_3
//RW-m3_h6 should use lag 7
gen d_rw_fa =e_RW_h6^2 - e_forecast_average_h6^2
ac d_rw_fa
corrgram  d_rw_fa
//RW-forecast averaging should use lag 7
gen d_rw_bg =e_RW_h6^2 - e_bate_granger_h6^2
ac d_rw_bg
corrgram  d_rw_bg
//RW-bates granger should use lag 7
gen d_rw_gr =e_RW_h6^2 - e_granger_raman_h6^2
ac d_rw_gr
corrgram  d_rw_gr
//RW-granger ramanathan should use lag 7
gen d_rw_bma =e_RW_h6^2 - e_bma_combined_h6^2
ac d_rw_bma
corrgram  d_rw_bma
//RW-bma should use lag 7
gen d_rw_waic =e_RW_h6^2 - e_waic_combined_h6^2
ac d_rw_waic
corrgram  d_rw_waic
//RW-waic should use lag 7


//dm test against AR1
gen d_ar1_ar12 =e_ar1_h6^2 - e_ar12_h6^2
ac d_ar1_ar12
corrgram  d_ar1_ar12
//AR1-AR12 should use lag 13
gen d_ar1_adl_1 =e_ar1_h6^2 - e_m2_h6^2
ac d_ar1_adl_1
corrgram  d_ar1_adl_1
//AR1-m2_h6 should use lag 13
gen d_ar1_adl_2 =e_ar1_h6^2 - e_m3m2_h6^2
ac d_ar1_adl_2
corrgram  d_ar1_adl_2
//AR1-m3m2_h6 should use lag 12
gen d_ar1_adl_3 =e_ar1_h6^2 - e_m3_h6^2
ac d_ar1_adl_3
corrgram  d_ar1_adl_3
//AR1-m3_h6 should use lag 12
gen d_ar1_fa =e_ar1_h6^2 - e_forecast_average_h6^2
ac d_ar1_fa
corrgram  d_ar1_fa
//AR1-forecast averaging should use lag 13
gen d_ar1_bg =e_ar1_h6^2 - e_bate_granger_h6^2
ac d_ar1_bg
corrgram  d_ar1_bg
//AR1-bates granger should use lag 13
gen d_ar1_gr =e_ar1_h6^2 - e_granger_raman_h6^2
ac d_ar1_gr
corrgram  d_ar1_gr
//AR1-granger ramanathan should use lag 13
gen d_ar1_bma =e_ar1_h6^2 - e_bma_combined_h6^2
ac d_ar1_bma
corrgram  d_ar1_bma
//AR1-bma should use lag 13
gen d_ar1_waic =e_ar1_h6^2 - e_waic_combined_h6^2
ac d_ar1_waic
corrgram  d_ar1_waic
//AR1-waic should use lag 13



//compiling results with dm test:
matrix rmse_matrix = J(11 , 3, 0) //11 rows with 3 columns

//compiling rmse (make sure its in order)
matrix rmse_matrix[1, 1] = rmse_RW_h6[1]
matrix rmse_matrix[2, 1] = rmse_ar1_h6[1]
matrix rmse_matrix[3, 1] = rmse_ar12_h6[1]
matrix rmse_matrix[4, 1] = rmse_m2_h6[1]
matrix rmse_matrix[5, 1] = rmse_m3m2_h6[1]
matrix rmse_matrix[6, 1] = rmse_m3_h6[1]
matrix rmse_matrix[7, 1] = rmse_forecast_average_h6[1]
matrix rmse_matrix[8, 1] = rmse_bate_granger_h6[1]
matrix rmse_matrix[9, 1] = rmse_granger_raman_h6[1]
matrix rmse_matrix[10, 1] = rmse_bma_combined_h6[1]
matrix rmse_matrix[11, 1] = rmse_waic_combined_h6[1]

//DM test against RW (make sure its in order):
qui dmariano true_values RW_h6 ar1_h6, crit(mse) kernel(bartlett) maxlag(1)
matrix rmse_matrix[2, 2] = r(p)
qui dmariano true_values RW_h6 ar12_h6, crit(mse) kernel(bartlett) maxlag(7)
matrix rmse_matrix[3, 2] = r(p)
qui dmariano true_values RW_h6 m2_h6, crit(mse) kernel(bartlett) maxlag(7)
matrix rmse_matrix[4, 2] = r(p)
qui dmariano true_values RW_h6 m3m2_h6, crit(mse) kernel(bartlett) maxlag(7)
matrix rmse_matrix[5, 2] = r(p)
qui dmariano true_values RW_h6 m3_h6, crit(mse) kernel(bartlett) maxlag(7)
matrix rmse_matrix[6, 2] = r(p)
qui dmariano true_values RW_h6 forecast_average_h6, crit(mse) kernel(bartlett) maxlag(7)
matrix rmse_matrix[7, 2] = r(p)
qui dmariano true_values RW_h6 bate_granger_h6, crit(mse) kernel(bartlett) maxlag(7)
matrix rmse_matrix[8, 2] = r(p)
qui dmariano true_values RW_h6 granger_raman_h6, crit(mse) kernel(bartlett) maxlag(7)
matrix rmse_matrix[9, 2] = r(p)
qui dmariano true_values RW_h6 bma_combined_h6, crit(mse) kernel(bartlett) maxlag(7)
matrix rmse_matrix[10, 2] = r(p)
qui dmariano true_values RW_h6 waic_combined_h6, crit(mse) kernel(bartlett) maxlag(7)
matrix rmse_matrix[11, 2] = r(p)

//DM test against AR(1) (make sure its in order):
qui dmariano true_values ar1_h6 ar12_h6, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[3, 3] = r(p)
qui dmariano true_values ar1_h6 m2_h6, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[4, 3] = r(p)
qui dmariano true_values ar1_h6 m3m2_h6, crit(mse) kernel(bartlett) maxlag(12)
matrix rmse_matrix[5, 3] = r(p)
qui dmariano true_values ar1_h6 m3_h6, crit(mse) kernel(bartlett) maxlag(12)
matrix rmse_matrix[6, 3] = r(p)
qui dmariano true_values ar1_h6 forecast_average_h6, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[7, 3] = r(p)
qui dmariano true_values ar1_h6 bate_granger_h6, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[8, 3] = r(p)
qui dmariano true_values ar1_h6 granger_raman_h6, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[9, 3] = r(p)
qui dmariano true_values ar1_h6 bma_combined_h6, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[10, 3] = r(p)
qui dmariano true_values ar1_h6 waic_combined_h6, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[11, 3] = r(p)

* Name the columns and rows for clarity 
matrix colnames rmse_matrix = "RMSE" "DM RW" "DM AR(1)"
matrix rownames rmse_matrix =  RW_h6 AR1_h6 AR12_h6 ADL1_h6 ADL2_h6 ADL3_h6 FA_h6 BG_h6 GR_h6 BMA_h6 WAIC_h6

matrix list rmse_matrix

putexcel set "6_step_results.xlsx", replace
putexcel A1=matrix(rmse_matrix), names










