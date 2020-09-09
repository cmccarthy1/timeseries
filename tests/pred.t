\l tm.q
.tm.loadfile`:init.q
\l tests/failMessage.q

\S 42
exogIntFuture   :1000 50#5000?1000
exogFloatFuture :1000 50#5000?1000f
exogMixedFuture :(1000 20#20000?1000),'(1000 20#20000?1000f),'(1000 10#10000?0b)

// Load files
fileList:`AR1`AR2`AR3`AR4`ARCH1`ARCH2`ARMA1`ARMA2`ARMA3`ARMA4`ARIMA1`ARIMA2,
         `ARIMA3`ARIMA4`SARIMA1`SARIMA2`SARIMA3`SARIMA4
loadFunc:{load hsym`$":tests/data/",x,string y}
loadFunc["fit/"]each fileList;
loadFunc["pred/pred"]each fileList;


// AR tests

.tm.AR.predict[AR1;();1000]~predAR1
.tm.AR.predict[AR2;exogFloatFuture;1000]~predAR2
.tm.AR.predict[AR3;exogIntFuture;1000]~predAR3
.tm.AR.predict[AR4;exogMixedFuture;1000]~predAR4

failingTest[.tm.AR.predict;(AR2;-1_'exogFloatFuture;1000);0b;"Test exog length does not match train exog length"]
failingTest[.tm.AR.predict;(AR3;-1_'exogIntFuture  ;1000);0b;"Test exog length does not match train exog length"]

// ARCH tests

.tm.ARCH.predict[ARCH1;1000]~predARCH1
.tm.ARCH.predict[ARCH2;1000]~predARCH2

// ARMA tests

.tm.ARMA.predict[ARMA1;();1000]~predARMA1
.tm.ARMA.predict[ARMA2;exogFloatFuture;1000]~predARMA2
.tm.ARMA.predict[ARMA3;exogIntFuture;1000]~predARMA3
.tm.ARMA.predict[ARMA4;exogMixedFuture;1000]~predARMA4

failingTest[.tm.ARMA.predict;(ARMA2;-1_'exogFloatFuture;1000);0b;"Test exog length does not match train exog length"                                                     ]
failingTest[.tm.ARMA.predict;(ARMA3;-1_'exogIntFuture  ;1000);0b;"Test exog length does not match train exog length"                                                     ]
failingTest[.tm.ARMA.predict;(AR1  ;()                 ;1000);0b;"The following required dictionary keys for 'mdl' are not provided: q_param, resid, estresid, pred_dict"]

// ARIMA tests

.tm.ARIMA.predict[ARIMA1;();1000]~predARIMA1
.tm.ARIMA.predict[ARIMA2;exogFloatFuture;1000]~predARIMA2
.tm.ARIMA.predict[ARIMA3;exogIntFuture;1000]~predARIMA3
.tm.ARIMA.predict[ARIMA4;exogMixedFuture;1000]~predARIMA4

failingTest[.tm.ARIMA.predict;(ARIMA2;-1_'exogFloatFuture;1000);0b;"Test exog length does not match train exog length"                       ]
failingTest[.tm.ARIMA.predict;(ARIMA3;-1_'exogIntFuture  ;1000);0b;"Test exog length does not match train exog length"                       ] 
failingTest[.tm.ARIMA.predict;(ARMA4 ;exogMixedFuture    ;1000);0b;"The following required dictionary keys for 'mdl' are not provided: origd"]

// SARIMA tests

.tm.SARIMA.predict[SARIMA1;();1000]~predSARIMA1
.tm.SARIMA.predict[SARIMA2;exogFloatFuture;1000]~predSARIMA2
.tm.SARIMA.predict[SARIMA3;exogIntFuture;1000]~predSARIMA3
.tm.SARIMA.predict[SARIMA4;exogMixedFuture;1000]~predSARIMA4

failingTest[.tm.SARIMA.predict;(SARIMA2;-1_'exogFloatFuture;1000);0b;"Test exog length does not match train exog length"                                         ]
failingTest[.tm.SARIMA.predict;(SARIMA3;-1_'exogIntFuture  ;1000);0b;"Test exog length does not match train exog length"                                         ]
failingTest[.tm.SARIMA.predict;(ARIMA2 ;exogFloatFuture    ;1000);0b;"The following required dictionary keys for 'mdl' are not provided: origs, P_param, Q_param"]


