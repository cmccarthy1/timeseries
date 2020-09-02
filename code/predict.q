\d .tm

// Prediction functionality for time-series models

// @kind function
// @category modelPredict
// @fileoverview Predictions based on an AutoRegressive model (AR)
// @param mdl  {dict} model parameters returned from fitting of an appropriate model
// @param exog {tab/num[][]/(::)} Exogenous variables, are additional variables which
//   required for application of model prediction 
// @param len  {integer} number of values to be predicted
// @return     {float[]} list of predicted values
ARpred:{[mdl;exog;len]
  exog:i.preddatacheck[mdl;exog];
  mdl[`pred_dict]:enlist[`p]!enlist count mdl`p_param;
  mdl[`estresid]:();
  mdl[`resid]:();
  i.predfunc[mdl;exog;len;i.sngpredAR]
  }

// @kind function
// @category modelPredict
// @fileoverview Predictions based on an AutoRegressive Moving Average model (ARMA)
// @param mdl  {dict} model parameters returned from fitting of an appropriate model
// @param exog {tab/num[][]/(::)} Exogenous variables, are additional variables which
//   required for application of model prediction 
// @param len  {integer} number of values to be predicted
// @return     {float[]} list of predicted values
// Predict future data using an ARMA model
/. r    > list of predicted values
ARMApred:{[mdl;exog;len]
  exog:i.preddatacheck[mdl;exog];
  i.predfunc[mdl;exog;len;i.sngpredARMA]
  }

// @kind function
// @category modelPredict
// @fileoverview Predictions based on an AutoRegressive Integrated Moving Average 
//   model (ARIMA)
// @param mdl  {dict} model parameters returned from fitting of an appropriate model
// @param exog {tab/num[][]/(::)} Exogenous variables, are additional variables which
//   required for application of model prediction 
// @param len  {integer} number of values to be predicted
// @return     {float[]} list of predicted values
ARIMApred:{[mdl;exog;len]
  exog:i.preddatacheck[mdl;exog];
  // Calculate predictions not accounting for differencing
  pred:i.predfunc[mdl;exog;len;i.sngpredARMA];
  dval:count mdl`origd;
  // Revert data to correct scale (remove differencing if previously applied)
  $[dval;dval _dval{sums x}/mdl[`origd],pred;pred]
  }

// @kind function
// @category modelPredict
// @fileoverview Predictions based on a Seasonal AutoRegressive Integrated Moving 
//   Average model (SARIMA)
// @param mdl  {dict} model parameters returned from fitting of an appropriate model
// @param exog {tab/num[][]/(::)} Exogenous variables, are additional variables which
//   required for application of model prediction 
// @param len  {integer} number of values to be predicted
// @return     {float[]} list of predicted values
SARIMApred:{[mdl;exog;len]
  exog:i.preddatacheck[mdl;exog];
  // Calculate predictions not accounting for differencing
  preds:$[count raze mdl[`pred_dict];
    i.predfunc[mdl;exog;len;i.sngpredSARMA];
    i.ARpred[mdl;exog;len]
    ];
  // Order of seasonal differencing originally applied
  sval:count mdl`origs;
  // if seasonal differenced, revert to original
  if[sval;preds:i.revseasdf[mdl[`origs];preds]];
  // Order of differencing originally applied
  dval:count mdl`origd;
  // Revert data to correct scale (remove differencing if previously applied)
  $[dval;dval _dval{sums x}/mdl[`origd],preds;preds]
  }


// @kind function
// @category modelPredict
// @fileoverview Predictions based on an AutoRegressive Conditional Heteroskedasticity 
//   model (ARCH)
// @param mdl  {dict} model parameters returned from fitting of an appropriate model
// @param len  {integer} number of values to be predicted
// @return     {float[]} list of predicted values
// Predict future volatility using an ARCH model
/. r    > list of predicted values
ARCHpred:{[mdl;len]
  // predict and return future values
  last{x>count y 1}[len;]i.sngpredARCH[mdl`params]/(mdl`resid;())
  }
