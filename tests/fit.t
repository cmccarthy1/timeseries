\l tm.q
.tm.loadfile`:init.q

// rounding function
round:{y*"j"$x%y}

// Time Series data
load`:tests/data/endogInt;
load`:tests/data/endogFloat;
load`:tests/data/exogInt;
load`:tests/data/exogFloat;
load`:tests/data/exogMixed;

// AR tests
keyAR:`params`tr_param`exog_param`p_param`lags

round[.tm.AR.fit[endogInt;();1;0b];.0001]~keyAR!(enlist -0.028;"f"$();"f"$();enlist -0.028;enlist 30f)
round[.tm.AR.fit[endogInt;exogFloat;3;1b];.0001]~keyAR!(49.5222 0.1873 -0.0182 -0.047 0.0152 0.066 -0.0697;enlist 49.5222;0.1873 -0.0182 -0.047;0.0152 0.066 -0.0697;50 37 30f)
round[.tm.AR.fit[endogFloat;exogInt;2;1b];.0001]~keyAR!(54.6993 0.0749 -0.1 0.1775 -0.0304 -0.1644;enlist 54.6993;0.0749 -0.1 0.1775;-0.0304 -0.1644;14.7755 27.4227)
round[.tm.AR.fit[endogFloat;exogMixed;4;0b];.0001]~keyAR!(0.0454 3.1558 -4.6593 0.2042 0.2432 0.1812 0.0239;"f"$();0.0454 3.1558 -4.6593;0.2042 0.2432 0.1812 0.0239;81.3248 12.9194 14.7755 27.4227)

// ARMA tests
keyARMA:`params`tr_param`exog_param`p_param`q_param`lags`resid`estresid`pred_dict

round[.tm.ARMA.fit[endogInt;();1;2;1b];.0001]~keyARMA!(50.3728 0.4404 0.1278 -0.5135;50.3728;"f"$();enlist 0.4404;0.1278 -0.5135;50 37 30f;29.7825 24.7519;0.0167 0.0968 -0.027;`p`q`tr!1 2 1f)
round[.tm.ARMA.fit[endogInt;exogFloat;2;1;0b];.0001]~keyARMA!(0.2822 0.0192 0.0121 0.1329 0.5791 -0.5331;0n;0.2822 0.0192 0.0121;0.1329 0.5791;enlist -0.5331;50 37 30f;enlist -15.5927;0.3365 0.0699 0.0785 0.1805 0.2249 0.1156;`p`q`tr!2 1 0f)
round[.tm.ARMA.fit[endogFloat;exogInt;1;1;0b];.0001]~keyARMA!(0.144 0.0353 0.2852 0.5637 -0.6212;0n;0.144 0.0353 0.2852;enlist 0.5637;enlist -0.6212;14.7755 27.4227;enlist -3.8496;0.265 0.0747 0.3909 0.1612 0.1005;`p`q`tr!1 1 0f)
round[.tm.ARMA.fit[endogFloat;exogMixed;3;2;1b];.0001]~keyARMA!(66.2214 -0.1215 1.1572 -3.4155 0.0127 -0.1823 -0.0864 0.1505 -0.086;66.2214;-0.1215 1.1572 -3.4155;0.0127 -0.1823 -0.0864;0.1505 -0.086;81.3248 12.9194 14.7755 27.4227;-49.1462 -7.58;0.0454 3.1558 -4.6593 0.2042 0.2432 0.1812 0.0239;`p`q`tr!3 2 1f)

//ARIMA tests
keyARIMA:keyARMA,`origd

round[.tm.ARIMA.fit[endogInt;();2;1;2;0b];.0001]~keyARIMA!(-0.3272 -1.6062 -0.8836 0.6788;0n;"f"$();-0.3272 -1.6062;-0.8836 0.6788; -34 -13 -7f; -14.9684 -24.3058;-0.1483 -0.4498 -0.8166;`p`q`tr!2 2 0f;enlist 30f)
round[.tm.ARIMA.fit[endogInt;exogFloat;1;1;1;1b];.0001]~keyARIMA!(4.3548 0.1768 -0.0847 -0.1627 -0.2495 -0.556;4.3548;0.1768 -0.0847 -0.1627;enlist -0.2495;enlist -0.556;-13 -7f;enlist -20.0261;0.2211 -0.0942 -0.1119 -0.305 -0.7431;`p`q`tr!1 1 1f;enlist 30f)
round[.tm.ARIMA.fit[endogFloat;exogInt;3;0;1;1b];.0001]~keyARIMA!(45.6925 0.0883 -0.0836 0.1464 -0.0015 -0.0387 0.0318 -0.1899;45.6925;0.0883 -0.0836 0.1464;-0.0015 -0.0387 0.0318;enlist -0.1899;81.3248 12.9194 14.7755 27.4227;enlist -14.4886;0.2221 0.0178 0.2803 0.1591 0.2112 0.0992 0.015;`p`q`tr!3 1 1f;"f"$())
round[.tm.ARIMA.fit[endogFloat;exogMixed;1;2;2;0b];.0001]~keyARIMA!(-0.0845 2.0728 -11.5966 -0.4124 -0.0614 -1.1562;0n;-0.0845 2.0728 -11.5966;enlist -0.4124;-0.0614 -1.1562;-118.5865 70.2615 10.7911;-27.2875 22.1198;0.0295 1.5475 -20.4809 -0.4777 -1.1034 -1.3508;`p`q`tr!1 2 0f;14.7755 -2.1282)

// SARIMA tests
keySARIMA:`params`tr_param`exog_param`p_param`q_param`P_param`Q_param`lags`resid`estresid`pred_dict`origd`origs
predDict:`p`q`P`Q`m`tr`seas_add_P`seas_add_Q`n
s1:`P`D`Q`m!1 0 2 5
s2:`P`D`Q`m!2 1 0 2
s3:`P`D`Q`m!2 1 1 3
s4:`P`D`Q`m!0 1 1 4

round[.tm.SARIMA.fit[endogInt;();1;1;1;0b;s1];.0001]~keySARIMA!(-0.1658 -0.6989 0.4511 -0.4817 -0.0476;0n;"f"$();enlist -0.1658;enlist -0.6989;enlist 0.4511;-0.4511 0.4817;-49 -2 58 -34 -13 -7f;10.2049 15.3563 -44.8139 -23.5985 35.745 -11.9463 -26.1234 39.9973 9.7958 -19.5731 -28.3947;-0.3361 -0.7667;predDict!(1f;1f;enlist 0f;0 5f;5f;0f;enlist 1f;1 6f;2f);enlist 30f;"f"$())
round[.tm.SARIMA.fit[endogInt;exogFloat;1;0;1;1b;s2];.0001]~keySARIMA!(3.3965 0.2321 -0.1337 -0.1636 -0.0285 -0.1174 -0.5352 -0.3;3.3965;0.2321 -0.1337 -0.1636;enlist -0.0285;enlist -0.1174;-0.5352 -0.3;"f"$();-51 56 24 -47 -20f;enlist -4.5175;0.2218 -0.1012 -0.1161 -0.4023 -0.0366;predDict!(1f;1f;0 2f;"f"$();2f;1f;1 3f;"f"$();2f);"f"$();37 30f)
round[.tm.SARIMA.fit[endogFloat;exogInt;1;2;0;0b;s3];.0001]~keySARIMA!(-0.1265 -0.1012 0.2456 -0.8023 0.0255 -0.1038 -0.2329;0n;-0.1265 -0.1012 0.2456;enlist -0.8023;"f"$();0.0255 -0.1038;enlist -0.0255;177.9025 -108.1708 23.6404 58.8881 -147.1471 157.5059 -94.739;-89.6168 22.892 41.9594;-0.1005 -0.2306 0.3378 -0.2708 -1.0234;predDict!(1f;0f;0 3f;enlist 0f;3f;0f;1 4f;"f"$();2f);14.7755 -2.1282;12.9194 14.7755 27.4227)
round[.tm.SARIMA.fit[endogFloat;exogMixed;2;1;1;0b;s4];.0001]~keySARIMA!(-0.077 0.6925 -0.1172 -0.0399 -0.1555 -0.7352 -0.609;0n;-0.077 0.6925 -0.1172;-0.0399 -0.1555;enlist -0.7352;"f"$();enlist 0.609;-71.7402 -30.0393 67.9962;-75.6087 37.1336 -43.439 -47.3541 6.4149; -0.0201 1.4339 -11.3271 0.0789 -0.3964 -0.786;predDict!(2f;1f;"f"$();enlist 0f;4f;0f;"f"$();enlist 1f;3f);enlist 27.4227;81.3248 12.9194 14.7755 27.4227)