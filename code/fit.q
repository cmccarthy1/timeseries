\d .tm

// Fitting functionality for time series models. 

// @kind function
// @category modelFit
// @fileoverview Fit an AutoRegressive model (AR)
// @param endog {num[]} Endogenous variable (time-series) from which to build a model
//   this is the target variable from which a value is to be predicted
// @param exog  {tab/num[][]/(::)} Exogenous variables, are additional variables which
//   may be accounted for to improve the model, if (::)/() this will be ignored
// @param lags  {integer} The number/order of time lags of the model
// @param trend {boolean} Is a trend line to be accounted for in fitting of model
// @return {dict} All information required to use a fit model for the prediction of
//   new values based on incoming data
ARfit:{[endog;exog;lags;trend]
  exog:i.fitdatacheck[endog;exog];
  // AR models require no non seasonal differencing steps
  i.differ[endog;0;()!()];
  // Estimate coefficients
  coeff:$[sum trend,count[exog];
    i.estparam[endog;exog;endog;`p`q`tr!lags,0,trend];
    i.durbin_lev[endog;lags]
    ];
  // Get lagged values needed for future predictions
  lagvals:neg[lags]#endog;
  // return dictionary with required info for predictions
  keyvals:`params`tr_param`exog_param`p_param`lags;
  params:(coeff;trend#coeff;coeff trend +til count exog 0;neg[lags]#coeff;lagvals);
  keyvals!params
  }

// @kind function
// @category modelFit
// @fileoverview Fit an AutoRegressive Moving Average model (ARMA)
// @param endog {num[]} Endogenous variable (time-series) from which to build a model
//   this is the target variable from which a value is to be predicted
// @param exog  {tab/num[][]/(::)} Exogenous variables, are additional variables which
//   may be accounted for to improve the model, if (::)/() this will be ignored
// @param lags  {integer} The number/order of time  lags of the model
// @param resid {integer} The number of residual errors to be accounted for
// @param trend {boolean} Is a trend line to be accounted for in fitting of model
// @return {dict} All information required to use a fit model for the prediction of
//   new values based on incoming data
ARMAfit:{[endog;exog;lags;resid;trend]
  $[resid~0;
    // if q = 0 then model is an AR model
    ARfit[endog;exog;lags;trend],`q_param`resid!(();());
    i.ARMAmodel[endog;exog;`p`q`tr!lags,resid,trend]]
  }

// @kind function
// @category modelFit
// @fileoverview Fit an AutoRegressive Integrated Moving Average model (ARIMA)
// @param endog {num[]} Endogenous variable (time-series) from which to build a model
//   this is the target variable from which a value is to be predicted
// @param exog  {tab/num[][]/(::)} Exogenous variables, are additional variables which
//   may be accounted for to improve the model, if (::)/() this will be ignored
// @param lags  {integer} The number/order of time  lags of the model
// @param diff  {integer} The order of time series differencing used in integration
// @param resid {integer} The number of residual errors to be accounted for
// @param trend {boolean} Is a trend line to be accounted for in fitting of model
// @return {dict} All information required to use a fit model for the prediction of
//   new values based on incoming data
ARIMAfit:{[endog;exog;lags;diff;resid;trend]
  exog:i.fitdatacheck[endog;exog];
  // Apply integration (non seasonal)
  I:i.differ[endog;diff;()!()];
  // Fit an ARMA model on the differenced time series
  mdl:ARMAfit[I;diff _exog;lags;resid;trend];
  // Retrieve the original data to be used when fitting on new data
  origData:neg[diff]#endog;
  // Produce the relevant differenced data for use in future predictions
  origDiff:enlist[`origd]!enlist diff{deltas x}/origData;
  // return relevant data
  mdl,origDiff
  }

// @kind function
// @category modelFit
// @fileoverview Fit a Seasonal AutoRegressive Integrated Moving Average model (SARIMA)
// @param endog {num[]} Endogenous variable (time-series) from which to build a model
//   this is the target variable from which a value is to be predicted
// @param exog  {tab/num[][]/(::)} Exogenous variables, are additional variables which
//   may be accounted for to improve the model, if (::)/() this will be ignored
// @param lags  {integer} The number/order of time  lags of the model
// @param diff  {integer} The order of time series differencing used in integration
// @param resid {integer} The number of residual errors to be accounted for
// @param trend {boolean} Is a trend line to be accounted for in fitting of model
// @param seas  {dict}    Is a dictionary containing required seasonal components
// @return {dict} All information required to use a fit model for the prediction of
//   new values based on incoming data
SARIMAfit:{[endog;exog;lags;diff;resid;trend;seas]
  i.dictCheck[seas;`P`Q`D`m;"seas"];
  // Apply error checking (exogenous data not converted to matrix?)
  exog:i.fitdatacheck[endog;exog];
  // Apply appropriate seasonal+non seasonal differencing
  I:i.differ[endog;diff;seas];
  // Create dictionary with p,q and seasonal components
  dict:`p`q`P`Q`m`tr!lags,resid,((1+til each seas[`P`Q])*seas[`m]),seas[`m],trend;
  // add additional seasonal components
  dict[`seas_add_P`seas_add_Q]:(raze'){1+til[x]+/:y}'[(lags;resid);dict`P`Q];
  // Generate data for regenerate data following differencing
  origDiffSeason:`origd`origs!(diff{deltas x}/neg[diff]#endog;neg[prd seas`D`m]#endog);
  // Apply SARMA model and postpend differenced original data
  i.SARMAmodel[I;exog;dict],origDiffSeason
  }

// @kind function
// @category modelFit
// @fileoverview Fit an AutoRegressive Condition Heteroscedasticity model (ARCH)
// @param endog {num[]} Endogenous variable (time-series) from which to build a model
//   this is the target variable from which a value is to be predicted
// @param exog  {tab/num[][]/(::)} Exogenous variables, are additional variables which
//   may be accounted for to improve the model, if (::)/() this will be ignored
// @param lags  {integer} The number/order of time  lags of the model
// @return {dict} All information required to use a fit model for the prediction of
//   new values based on incoming data
// Fit an ARCH model
/. r    > the model parameters and data needed for future predictions
ARCHfit:{[endog;exog;lags]
  exog:i.fitdatacheck[endog;exog];
  // Retrieve squared error
  sqer:err*err:i.esterrs[endog;exog;lags]`err;
  // Using the resid errorrs calculate coefficients
  coeff:i.estparam[sqer;();sqer;`p`q`tr!lags,0,1b];
  // Get lagged values needed for future predictions
  resid:neg[lags]#sqer;
  // return dictionary with required info for predictions
  keyvals:`params`tr_param`p_param`resid;
  params:(coeff;coeff[0];1_coeff;resid);
  keyvals!params
  }
