\l tm.q
.tm.loadfile`:init.q
\l tests/failMessage.q

\S 42
endogInt  :10000?1000
endogFloat:10000?1000f
exogInt   :10000 50#50000?1000
exogFloat :10000 50#50000?1000f
exogMixed :(10000 20#20000?1000),'(10000 20#20000?1000f),'(10000 10#10000?0b)
residInt  :10000?1000
residFloat:10000?1000f

// Load files
fileList:`AR1`AR2`AR3`AR4`ARCH1`ARCH2`ARMA1`ARMA2`ARMA3`ARMA4`ARIMA1`ARIMA2,
         `ARIMA3`ARIMA4`SARIMA1`SARIMA2`SARIMA3`SARIMA4`nonStat
{load hsym`$":tests/data/fit/",string x}each fileList;

// AR tests
.tm.AR.fit[endogInt  ;()       ;1;0b]~AR1
.tm.AR.fit[endogInt  ;exogFloat;3;1b]~AR2
.tm.AR.fit[endogFloat;exogInt  ;2;1b]~AR3
.tm.AR.fit[endogFloat;exogMixed;4;0b]~AR4

failingTest[.tm.AR.fit;(endogInt  ;5000#exogInt  ;1;1b);0b;"Endog length less than length"]
failingTest[.tm.AR.fit;(endogFloat;5000#exogFloat;1;1b);0b;"Endog length less than length"]


// ARMA tests
.tm.ARMA.fit[endogInt  ;()       ;1;2;1b]~ARMA1
.tm.ARMA.fit[endogInt  ;exogFloat;2;1;0b]~ARMA2
.tm.ARMA.fit[endogFloat;exogInt  ;1;1;0b]~ARMA3
.tm.ARMA.fit[endogFloat;exogMixed;3;2;1b]~ARMA4

failingTest[.tm.ARMA.fit;(endogInt  ;5000#exogInt  ;2;1;0b);0b;"Endog length less than length"                     ]
failingTest[.tm.ARMA.fit;(endogFloat;5000#exogFloat;2;1;0b);0b;"Endog length less than length"                     ]


// ARCH tests
.tm.ARCH.fit[residInt  ;3]~ARCH1
.tm.ARCH.fit[residFloat;1]~ARCH2


// ARIMA tests
.tm.ARIMA.fit[endogInt  ;()       ;2;1;2;0b]~ARIMA1
.tm.ARIMA.fit[endogInt  ;exogFloat;1;1;1;1b]~ARIMA2
.tm.ARIMA.fit[endogFloat;exogInt  ;3;0;1;1b]~ARIMA3
.tm.ARIMA.fit[endogFloat;exogMixed;1;2;2;0b]~ARIMA4

failingTest[.tm.ARIMA.fit;(endogInt  ;5000#exogInt  ;1;1;1;1b);0b;"Endog length less than length"                     ]
failingTest[.tm.ARIMA.fit;(endogFloat;5000#exogFloat;1;1;1;1b);0b;"Endog length less than length"                     ]
failingTest[.tm.ARIMA.fit;(nonStat   ;()            ;1;0;1;1b);0b;"Time series not stationary, try another value of d"]

// SARIMA tests
s1:`P`D`Q`m!1 0 2 5
s2:`P`D`Q`m!2 1 0 10
s3:`P`D`Q`m!2 1 1 30
s4:`P`D`Q`m!0 1 1 20

.tm.SARIMA.fit[endogInt  ;()       ;1;1;1;0b;s1]~SARIMA1
.tm.SARIMA.fit[endogInt  ;exogFloat;1;0;1;1b;s2]~SARIMA2
.tm.SARIMA.fit[endogFloat;exogInt  ;1;2;0;0b;s3]~SARIMA3
.tm.SARIMA.fit[endogFloat;exogMixed;2;1;1;0b;s4]~SARIMA4

failingTest[.tm.SARIMA.fit;(endogInt  ;5000#exogInt  ;2;0;1;1b;s1);0b;"Endog length less than length"                     ]
failingTest[.tm.SARIMA.fit;(endogFloat;5000#exogFloat;2;0;1;1b;s1);0b;"Endog length less than length"                     ]
failingTest[.tm.SARIMA.fit;(nonStat   ;()            ;2;0;0;1b;s1);0b;"Time series not stationary, try another value of d"]

