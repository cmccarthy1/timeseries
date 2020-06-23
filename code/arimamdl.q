// Load relevant utilities

\d .tm

// Fitting functionality for AR/ARMA/ARIMA models

/* endog = Endogenous variable (time-series)
/* exog  = Exogenous variables, if (::)/() then exogenous variables ignored
/* p     = Number of lags
/* q     = Number of residual errors
/* d     = Order of differencing being used
/* tr    = Is a trend to be accounted for

// Fit a basic AutoRegressive (AR) model
/. r the model parameters and data needed for future predictions
ARfit:{[endog;exog;p;tr]
  // check that exogenous variables are appropriate
  if[exog~(::);exog:()];
  if[not[()~exog]&count[endog]<>count exog;i.err.len[]];
  if[98h~type exog;exog:"f"$i.mat exog];
  // Check that the dataset must be stationary
  if[not i.stat[endog];i.err.stat[]];
  // Estimate coefficients
  coeff:i.estparam[endog;exog;endog;p;0;tr];
  // Get lagged values needed for future predictions
  lagvals:neg[p]#endog;
  // return dictionary with required info for predictions
  keyvals:`params`tr_param`exog_param`p_param`lags;
  params:(coeff;tr#coeff;coeff tr _til count exog;neg[p]#coeff;lagvals);
  keyvals!params
  }

// Fit an AutoRegressive Moving Average (ARMA) model
/. r > the model parameters and relevant data to accommodate future predictions
ARMAfit:{[endog;exog;p;q;tr]
  $[q~0;
    // if q = 0 then model is and AR model
    ARfit[endog;exog;p;tr],`q_param`resid!(();());
    i.ARMAmdl[endog;exog;p;q;tr]]
  }

// Fit an AutoRegressive Integrated Moving Average (ARIMA) model
/. r > the model parameters and relevant data to accommodate future predictions
ARIMAfit:{[endog;exog;p;d;q;tr]
  // Check that exogenous variables are appropriate
  if[exog~(::);exog:()];
  if[not[()~exog]&count[endog]<>count exog;i.err.len[]];
  // 'Integrate' the time series through differencing
  I:i.diff[endog;d];
  // Check the dataset is stationary
  if[not i.stat I;i.err.stat[]];
  // Fit an ARMA model on the differenced time series
  mdl:ARMAfit[I;exog;p;q;tr];
  // Produce the relevant differenced data for use in future predictions
  orig_diff:enlist[`origd]!enlist d{deltas x}/neg[d]#endog;
  // return relevant data
  mdl,orig_diff
  }


// Prediction functionality for AR/ARMA/ARIMA models

// Predict future data using an AR model
/* mdl  = model parameters returned from fitted AR model
/* exog = any exogenous values that are required for the fitting (endogenous values provided in model)
/* len  = number of values that are to be predicted
/. r    > list of predicted values
ARpred:{[mdl;exog;len]
  // allow null to be provided as exogenous variable
  if[exog~(::);exog:()];
  // convert exogenous variable to a matrix if required
  if[98h~type exog;exog:"f"$i.mat exog];
  // predict and return future values
  last{x>count y 2}[len;]i.sngpred[mdl`params;exog;0;;()]/(mdl`lags;();())
  }

ARMApred:{[mdl;exog;len]
  // allow null to be provided as exogenous variable
  if[exog~(::);exog:()];
  // convert exogenous variable to a matrix if required
  if[98h~type exog;exog:"f"$i.mat exog];
  // Calculate the difference between number of predicted values (mdl`lags) and number of coefficients
  pval:count[mdl`lag]-count mdl`p_param;
  // if any residual values are present then use these in prediction otherwise use ARpred
  $[count mdl`resid;
    last{x>count y 2}[len;]i.sngpred[mdl`params;exog;pval;;mdl`estresid]/(mdl`lags;mdl`resid;());
    ARpred[mdl;exog;len]]
  }

ARIMApred:{[mdl;exog;len]
  // Predict values
  pred:ARMApred[mdl;exog;len];
  // Order of differencing originally applied
  dval:count mdl`origd;
  // Revert data to correct scale (remove differencing if previously applied)
  $[dval;dval _dval{count[x] msum x}/mdl[`origd],::;]preds
  }


// Summary of the stationarity of each column of a multivariate time series or a single vector
/* tseries = table containing the time series columns in question or a vector of values comprising
/*           a time series of interest, note the entries should in each case be numeric data types
/. r > a keyed table with some of the most informative outputs from the python adfuller test

stationary:{[dset]
  dtype:type dset;
  // Names to be provided to form the key for the return table
  keynames:$[99h=dtype;key dset;98h=dtype;cols dset;enlist`data];
  // Column names associated with the returns from the augmented dickey fuller test
  dcols:`ADFstat`pvalue`stationary,`$raze each"CriticalValue_",/:string(1;5;10),\:"%";
  scores:i.statscores[dset;dtype];
  keynames!flip dcols!scores
  }

