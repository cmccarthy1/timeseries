\l tm.q
.tm.loadfile`:init.q

// rounding function
round:{y*"j"$x%y}

// Time Series data
\l tm.q
.tm.loadfile`:init.q

// Time Series data
load`:tests/data/endogInt;
load`:tests/data/endogIntFuture;
load`:tests/data/endogFloat;
load`:tests/data/endogFloatFuture;
load`:tests/data/exogInt;
load`:tests/data/exogIntFuture;
load`:tests/data/exogFloat;
load`:tests/data/exogFloatFuture;
load`:tests/data/exogMixed;
load`:tests/data/exogMixedFuture;

// Stationality

// Staionary columns
dcols:`ADFstat`pvalue`stationary,`$raze each"CriticalValue_",/:string(1;5;10),\:"%";

// Test return of stationality functions
cols[value .tm.stationality[endogInt]]~dcols
cols[value .tm.stationality[endogFloat]]~dcols

// aicparam 

// Set up parameters
dictKeys :`endog`exog
paramKeys:`p`d`q`tr
aicKeys  :`p`d`q`tr`score

trainDict1:dictKeys!(endogInt  ;()       )
trainDict2:dictKeys!(endogInt  ;exogFloat)
trainDict3:dictKeys!(endogFloat;exogInt  )
trainDict4:dictKeys!(endogFloat;exogMixed)

testDict1:dictKeys!(endogIntFuture  ;()             )
testDict2:dictKeys!(endogIntFuture  ;exogFloatFuture)
testDict3:dictKeys!(endogFloatFuture;exogIntFuture  )
testDict4:dictKeys!(endogFloatFuture;exogMixedFuture)

params:paramKeys!(1 1 3 2;1 0 1 0;1 0 2 0;0011b)

// Test return of aicparam
round[.tm.aicParam[trainDict1;testDict1;10;params];.001]~aicKeys!1 1 1 0 77.013
round[.tm.aicParam[trainDict2;testDict2;10;params];.001]~aicKeys!1 0 0 0 70.563
round[.tm.aicParam[trainDict3;testDict3;10;params];.001]~aicKeys!1 0 0 0 74.019
round[.tm.aicParam[trainDict4;testDict4;10;params];.001]~aicKeys!1 0 0 0 74.024


// Feature Exraction time Series tables

// Set up tables
ts_tab:([]"p"$"d"$til 10;10?10f;10?100)

// Windowed features
cols[.tm.tsWindow[ts_tab;`x1`x2;`max`min`avg;2 3]]~`x`x1`x2`max_2_x1`max_2_x2`max_3_x1`max_3_x2`min_2_x1`min_2_x2`min_3_x1`min_3_x2`avg_2_x1`avg_2_x2`avg_3_x1`avg_3_x2
count[.tm.tsWindow[ts_tab;`x1`x2;enlist`max;2 4]]~6

// Lagged Features
cols[.tm.tsLag[ts_tab;`x1`x2;enlist 3]]~`x`x1`x2`x1_xprev_3`x2_xprev_3
count[.tm.tsLag[ts_tab;enlist`x2;enlist 3]]~10