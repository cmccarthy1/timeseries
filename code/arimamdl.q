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
/. r > the model parameters and data needed for future predictions
ARfit:{[endog;exog;p;tr]
  // check that exogenous variables are appropriate
  if[exog~(::);exog:()];
  if[not[()~exog]&(count[endog])>count exog;i.err.len[]];
  if[98h~type exog;exog:"f"$i.mat exog];
  // Check that the dataset must be stationary
  if[not i.stat[endog];i.err.stat[]];
  // Estimate coefficients
  coeff:i.estparam[endog;exog;endog;`p`q!p,0;tr];
  // Get lagged values needed for future predictions
  lagvals:neg[p]#endog;
  // return dictionary with required info for predictions
  keyvals:`params`tr_param`exog_param`p_param`lags;
  params:(coeff;tr#coeff;coeff tr _til count exog 0;neg[p]#coeff;lagvals);
  keyvals!params
  }

// Fit an AutoRegressive Moving Average (ARMA) model
/. r    > the model parameters and data needed for future predictions
ARMAfit:{[endog;exog;p;q;tr]
  $[q~0;
    // if q = 0 then model is and AR model
    ARfit[endog;exog;p;tr],`q_param`resid!(();());
    i.SARMAmdl[endog;exog;`p`q!p,q;tr;"ARMA"]]
  }

// Fit an AutoRegressive Integrated Moving Average (ARIMA) model
/. r    > the model parameters and data needed for future predictions
ARIMAfit:{[endog;exog;p;d;q;tr]
  // Check that exogenous variables are appropriate
  if[exog~(::);exog:()];
  if[not[()~exog]&count[endog]>count exog;i.err.len[]];
  // 'Integrate' the time series through differencing
  I:i.diff[endog;d];
  // Check the dataset is stationary
  if[not i.stat I;i.err.stat[]];
  // Fit an ARMA model on the differenced time series
  mdl:ARMAfit[I;d _exog;p;q;tr];
  // Produce the relevant differenced data for use in future predictions
  orig_diff:enlist[`origd]!enlist d{deltas x}/neg[d]#endog;
  // return relevant data
  mdl,orig_diff
  }

// Fit a Seasonal ARIMA (SARIMA) model
/*s      = a dictionary of seasonal components
/. r     > the model parameters and data needed for future predictions
SARIMAfit:{[endog;exog;p;d;q;tr;s]
 // chk length
 if[exog~(::);exog:()];
 if[not[()~exog]&count[endog]>count exog;i.err.len[]];
 // diff model
 I:i.diff[endog;d];
 // Seasonality Diff
 if[s`D;I:s[`D]i.sdiff[s`m]/I];
 // do stationary check
 if[not i.stat[I];i.err.stat];
 // Create dictionary with p,q and seasonal components
 dict:`p`q`P`Q!p,q,(1+til each s[`P`Q])*s[`m];
 // run ARMA model
 i.SARMAmdl[I;exog;dict;tr;"SAR"],`origd`origs!(d{deltas x}/neg[d] #endog;neg[s[`D]*s`m]#endog)}

// Fit an ARCH model
/. r    > the model parameters and data needed for future predictions
ARCHfit:{[endog;exog;p]
 // check that exogenous variables are appropriate 
 if[not[()~exog]&(count[endog])>count exog;i.err.len[]];
 // convert exon table to matrix 
 if[98h~type exog;exog:"f"$i.mat[exog]];
 errs:i.esterrs[endog;exog;p];
 sqer:errs*errs:errs`err;
 // Using the resid errorrs calculate coeffi
 coeff:i.estparam[sqer;();sqer;`p`q!p,0;1b];
 // Get lagged values needed for future predictions
 resid:neg[p]#sqer;
 // return dictionary with required info for predictions
 keyvals:`params`tr_param`p_param`resid;
 params:(coeff;coeff[0];1_coeff;resid);
 keyvals!params
 }

// Prediction functionality for AR/ARMA/ARIMA/SARIMA models
/* mdl  = model parameters returned from fitted AR model
/* exog = any exogenous values that are required for the fitting (endogenous values provided in model)
/* len  = number of values that are to be predicted

// Predict future data using an AR model
/. r    > list of predicted values
ARpred:{[mdl;exog;len]
  // allow null to be provided as exogenous variable
  if[exog~(::);exog:()];
  // convert exogenous variable to a matrix if required
  if[98h~type exog;exog:"f"$i.mat exog];
  // predict and return future values
  last{x>count y 2}[len;]i.sngpred[mdl`params;exog;enlist[`p]!enlist count mdl`p_param;;();"ARMA"]/(mdl`lags;();())
  }

// Predict future data using an ARMA model
/. r    > list of predicted values
ARMApred:{[mdl;exog;len]
  i.SARMApred[mdl;exog;len;"ARMA"]
  }

// Predict future data using an ARIMA model
/. r    > list of predicted values
ARIMApred:{[mdl;exog;len]
  // Predict values
  pred:i.SARMApred[mdl;exog;len;"ARMA"];
  // Order of differencing originally applied
  dval:count mdl`origd;
  // Revert data to correct scale (remove differencing if previously applied)
  $[dval;dval _dval{sums x}/mdl[`origd],pred;pred]
  }

// Predict future data using a SARIMA model
/. r    > list of predicted values
SARIMApred:{[mdl;exog;len]
 // predict values
 if[98h~type exog;exog:"f"$.tm.i.mat[exog]];
 // if MA=0 then use ARpred
 preds:i.SARMApred[mdl;exog;len;"SAR"];
 / Order of seasonal differencing originally applied
 sval:count mdl`origs;
 // if seasonal differenced, revert to original
 if[sval;preds:i.revseasdf[mdl[`origs];preds]];
 // Order of differencing originally applied
  dval:count mdl`origd;
 // Revert data to correct scale (remove differencing if previously applied)
 $[dval;dval _dval{sums x}/mdl[`origd],preds;preds]}

// Predict future volatility using an ARCH model
/. r    > list of predicted values
ARCHpred:{[mdl;len]
  // predict and return future values
  last{x>count y 1}[len;]{[params;d]
  preds:params[0]+d[0] mmu 1_params;
  ((1_d[0]),preds;d[1],preds)
  }[mdl`params]/(mdl`resid;())
  }

// Summary of the stationarity of each column of a multivariate time series or a single vector
/* dset    = table containing the time series columns in question or a vector of values comprising
/*           a time series of interest, note the entries should in each case be numeric data types
/. r       > a keyed table with some of the most informative outputs from the python adfuller test
stationary:{[dset]
  dtype:type dset;
  // Names to be provided to form the key for the return table
  keynames:$[99h=dtype;key dset;98h=dtype;cols dset;enlist`data];
  // Column names associated with the returns from the augmented dickey fuller test
  dcols:`ADFstat`pvalue`stationary,`$raze each"CriticalValue_",/:string(1;5;10),\:"%";
  scores:i.statscores[dset;dtype];
  keynames!flip dcols!scores
  }

/Return best parameters based on lowes AIC score
/* train  = training data dictionary with keys `endog`exog
/* test   = testing data dictionry  with keys `engod`exog
/* len    = how many steps ahead to predict
/* dict   = dictionary of parameters to try
/. r      > params that scored the best
aicparam:{[train;test;len;dict]
 // get aic scores for each set of params
 scores:i.aicfitscore[train;test;len;]each flip dict;
 // return best value
 dict[;scores?bsc],enlist[`score]!enlist bsc:min scores}

