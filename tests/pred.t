\l tm.q
.tm.loadfile`:init.q

// rounding function
round:{y*"j"$x%y}

// Time Series data
load`:tests/data/endogInt;
load`:tests/data/endogFloat;
load`:tests/data/exogInt;
load`:tests/data/exogIntFuture;
load`:tests/data/exogFloat;
load`:tests/data/exogFloatFuture;
load`:tests/data/exogMixed;
load`:tests/data/exogMixedFuture;

// AR tests

// Fit AR models using data
ARmdl1:.tm.AR.fit[endogInt  ;()       ;1;0b]
ARmdl2:.tm.AR.fit[endogInt  ;exogFloat;3;1b]
ARmdl3:.tm.AR.fit[endogFloat;exogInt  ;2;1b]
ARmdl4:.tm.AR.fit[endogFloat;exogMixed;2;0b]


// Test predict function
round[.tm.AR.predict[ARmdl1;();10];.0001]~-0.8413 0.0236 -0.0007 0 0 0 0 0 0 0
round[.tm.AR.predict[ARmdl2;exogFloatFuture;10];.0001]~64.5429 50.5236 53.6283 58.4118 47.2557 57.2833 56.2306 52.3903 55.6782 57.8418
round[.tm.AR.predict[ARmdl3;exogIntFuture;10];.0001]~62.5472 45.598 44.9287 54.9193 48.415 48.8699 43.213 52.4084 54.0587 55.7673
round[.tm.AR.predict[ARmdl4;exogMixedFuture;10];.0001]~52.2428 27.9727 38.9342 41.4862 31.9208 34.2384 52.3429 32.29 53.1631 36.8018

// ARMA tests

// Fit ARMA models using data
ARMAmdl1:.tm.ARMA.fit[endogInt ;()        ;1;2;0b]
ARMAmdl2:.tm.ARMA.fit[endogInt ;exogFloat ;3;1;1b]
ARMAmdl3:.tm.ARMA.fit[endogFloat;exogInt  ;2;2;1b]
ARMAmdl4:.tm.ARMA.fit[endogFloat;exogMixed;2;1;0b]

// Test predict function
round[.tm.ARMA.predict[ARMAmdl1;();10];.0001]~37.9415 31.3854 31.109 31.1962 28.9951 28.4025 27.4961 26.3438 25.5458 24.6419
round[.tm.ARMA.predict[ARMAmdl2;exogFloatFuture;10];.0001]~61.5181 52.4786 50.7939 56.2504 44.9437 54.5128 54.9402 50.0986 52.559 56.7866
round[.tm.ARMA.predict[ARMAmdl3;exogIntFuture;10];.0001]~60.9679 46.2905 35.0551 54.2071 58.9245 45.7144 44.2897 46.2296 57.3992 53.1521
round[.tm.ARMA.predict[ARMAmdl4;exogMixedFuture;10];.0001]~44.6786 29.4319 25.6966 41.6506 34.4415 26.2086 41.6739 44.8297 41.8547 41.405

// ARIMA tests

// Fit ARIMA models using data
ARIMAmdl1:.tm.ARIMA.fit[endogInt ;()        ;2;1;1;0b]
ARIMAmdl2:.tm.ARIMA.fit[endogInt ;exogFloat ;3;1;0;1b]
ARIMAmdl3:.tm.ARIMA.fit[endogFloat;exogInt  ;2;2;1;1b]
ARIMAmdl4:.tm.ARIMA.fit[endogFloat;exogMixed;1;0;1;0b]

// Test predict function
round[.tm.ARIMA.predict[ARIMAmdl1;();10];.0001]~47.1647 42.6365 40.4593 40.7045 41.3495 41.665 41.1738 41.2324 41.3001 41.3161
round[.tm.ARIMA.predict[ARIMAmdl2;exogFloatFuture;10];.0001]~49.5984 45.0963 44.4946 49.6426 38.4557 46.8803 45.6338 38.6672 38.7323 43.7127
round[.tm.ARIMA.predict[ARIMAmdl3;exogIntFuture;10];.0001]~-10.7876   -8.0443  -44.4119  -37.7519  -65.1124  -93.7279 -116.76 -140.7334 -157.2691 -180.473
round[.tm.ARIMA.predict[ARIMAmdl4;exogMixedFuture;10];.0001]~36.0012 39.8854 31.2191 44.2636 35.8543 32.1093 45.9513 48.7715 44.0065 46.0751

// SARIMA tests

// Seasonal Dictionaries
s1:`P`D`Q`m!1 0 2 5
s2:`P`D`Q`m!2 1 0 2
s3:`P`D`Q`m!2 1 1 3
s4:`P`D`Q`m!0 1 1 4

// Fit SARIMA models using data
SARIMAmdl1:.tm.SARIMA.fit[endogInt ;()        ;1;1;1;0b;s1]
SARIMAmdl2:.tm.SARIMA.fit[endogInt ;exogFloat ;2;0;1;1b;s2]
SARIMAmdl3:.tm.SARIMA.fit[endogFloat;exogInt  ;1;0;0;1b;s3]
SARIMAmdl4:.tm.SARIMA.fit[endogFloat;exogMixed;3;2;2;0b;s4]

// Test predict function
round[.tm.SARIMA.predict[SARIMAmdl1;();10];.0001]~29.0054 11.5174 69.5879 36.8991 11.5641 32.2171 62.3462 5.4308 49.0707 20.4636
round[.tm.SARIMA.predict[SARIMAmdl2;exogFloatFuture;10];.0001]~-3.7282  111.3508  -73.8138  152.1897  -93.3454  172.2161  -38.0381 52.5464   84.3373 -108.2345
round[.tm.SARIMA.predict[SARIMAmdl3;exogIntFuture;10];.0001]~-37.6009  32.1923  27.1999 -66.9061 124.6995  29.1426 -33.8472 211.8595 -12.8421 -19.6088
round[.tm.SARIMA.predict[SARIMAmdl4;exogMixedFuture;10];.0001]~-6.1482  65.6318  67.3657 184.1249 208.5847 319.2122 369.1547 555.6313 676.709  808.0481