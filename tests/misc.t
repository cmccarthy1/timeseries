\l tm.q
.tm.loadfile`:init.q

\S 42

// Training data
endogInt  :10000?1000
endogFloat:10000?1000f
exogInt   :10000 50#50000?1000
exogFloat :10000 50#50000?1000f
exogMixed :(10000 20#20000?1000),'(10000 20#20000?1000f),'(10000 10#10000?0b)

// Testing data
endogIntFuture  :1000?1000
endogFloatFuture:1000?1000f
exogIntFuture   :1000 50#5000?1000
exogFloatFuture :1000 50#5000?1000f
exogMixedFuture :(1000 20#20000?1000),'(1000 20#20000?1000f),'(1000 10#10000?0b)


// Load files
fileList:`stationalityTab1`stationalityTab2`aicScore1`aicScore2`aicScore3`aicScore4,
         `windowTab1`windowTab2`lagTab1`lagTab2
{load hsym`$":tests/data/misc/",string x}each fileList;


// Stationality

.tm.stationality[endogInt  ]~stationalityTab1
.tm.stationality[endogFloat]~stationalityTab2

// aicparam 

// Set up parameters
dictKeys :`endog`exog
paramKeys:`p`d`q`tr

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
.tm.aicParam[trainDict1;testDict1;1000;params]~aicScore1
.tm.aicParam[trainDict2;testDict2;1000;params]~aicScore2
.tm.aicParam[trainDict3;testDict3;1000;params]~aicScore3
.tm.aicParam[trainDict4;testDict4;1000;params]~aicScore4

// Feature Exraction time Series tables

// Set up tables
ts_tab:([]"p"$"d"$til 1000;1000?10f;1000?100;1000?1f;1000?1000)

// Test windowed features
.tm.tsWindow[ts_tab;`x1`x2`x3`x4;`max`min`avg;2 3     ]~windowTab1
.tm.tsWindow[ts_tab;`x1`x2      ;`min`avg    ;enlist 2]~windowTab2

// Test lagged features
.tm.tsLag[ts_tab;`x1`x2`x4;enlist 3]~lagTab1
.tm.tsLag[ts_tab;enlist`x1;2 4 6   ]~lagTab2
