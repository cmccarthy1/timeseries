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
AR.predict:{[mdl;exog;len]
  i.dictCheck[mdl;i.AR.keyList;"mdl"];
  exog:i.predDataCheck[mdl;exog];
  mdl[`pred_dict]:enlist[`p]!enlist count mdl`p_param;
  mdl[`estresid]:();
  mdl[`resid]:();
  i.predictFunction[mdl;exog;len;i.AR.singlePredict]
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
ARMA.predict:{[mdl;exog;len]
  i.dictCheck[mdl;i.ARMA.keyList;"mdl"];
  exog:i.predDataCheck[mdl;exog];
  i.predictFunction[mdl;exog;len;i.ARMA.singlePredict]
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
ARIMA.predict:{[mdl;exog;len]
  i.dictCheck[mdl;i.ARIMA.keyList;"mdl"];
  exog:i.predDataCheck[mdl;exog];
  // Calculate predictions not accounting for differencing
  pred:i.predictFunction[mdl;exog;len;i.ARMA.singlePredict];
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
SARIMA.predict:{[mdl;exog;len]
  i.dictCheck[mdl;i.SARIMA.keyList;"mdl"];
  exog:i.predDataCheck[mdl;exog];
  // Calculate predictions not accounting for differencing
  preds:$[count raze mdl[`pred_dict];
    i.predictFunction[mdl;exog;len;i.SARMA.singlePredict];
    i.AR.predict[mdl;exog;len]
    ];
  // Order of seasonal differencing originally applied
  sval:count mdl`origs;
  // if seasonal differenced, revert to original
  if[sval;preds:i.reverseSeasonDiff[mdl[`origs];preds]];
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
ARCH.predict:{[mdl;len]
  i.dictCheck[mdl;i.ARCH.keyList;"mdl"];
  // predict and return future values
  last{x>count y 1}[len;]i.ARCH.singlePredict[mdl`params]/(mdl`resid;())
  }
