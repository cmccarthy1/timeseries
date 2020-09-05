\l ml/ml.q
.ml.loadfile`:init.q

// Namespace appropriately
\d .tm

// optimization function based off scipy https://docs.scipy.org/doc/scipy/reference/optimize.minimize-bfgs.html#optimize-minimize-bfgs
// Derivation of Formula http://www.bioinfo.org.cn/~wangchao/maa/Numerical_Optimization.pdf

// epsilon value for gradient fnc hardcoded for the moment,
// potentially should be changable by the user in the future
eps:1.49e-8;

funcEval:{[f;x0;args]
  $[any args~/:((::);());f x0;
    99h=type args;f[x0]. value args;
    '"args must be a dictionary, (::) or ()"
   ]
  }

updDefault:{[params]
  returnKeys:`norm`maxiter`gtol`geps`display`stepsize`c1`c2;
  returnVals:(0W;0W;1e-4;1.49e-8;0b;0W;1e-4;0.9);
  returnDict:returnKeys!returnVals;
  if[99h<>type params;params:()!()];
  returnDict,params
  }

vecNorm:{[vec;ord]
  if[-7h<>type ord;'"ord must be +/- infinity or a long atom"];
  $[ 0W~ord;max abs vec;
    -0W~ord;min abs vec;
    sum[abs[vec]xexp ord]xexp 1%ord
  ]
  }

stopCondition:{[dict;params]
  (dict[`fk] < dict`prev_fk) & (not any null dict`xk) & 
  (params[`maxiter] > dict`idx ) & (params[`gtol] < dict`gnorm)
  }


i.gradEval:{[fk;func;xk;args;eps;idx]
  if[(::)~fk;fk:funcEval[func;xk;args]];
  // increment function optimisation values by epsilon
  xk[idx]+:eps;
  // Evaluate the gradient
  (funcEval[func;xk;args]-fk)%eps
  }

// function with step value applied
/* f     = function to be minimized
/* pk    = step direction
/* alpha = step applied
/* d     > dictionary of x variables and argument
/. r - returns value of function with step value applied
i.phi:{[f;pk;alpha;xk;args]
  xk+:alpha*pk;
  funcEval[f;xk;args]
  }

// derivative of function with step value included
/* f     = function to be minimized
/* pk    = step direction
/* alpha = step applied
/* d     = dictionary of x variables and argument
/. r     > returns dictionary of gradiant and scalar derivative
i.derphi:{[func;eps;pk;alpha;xk;args;fk]
  // increment xk by a small step size
  xk+:alpha*pk;
  // get gradient at the new position
  gval:i.grad[func;xk;args;fk;eps];
  derval:gval mmu pk;
  `grad`derval!(gval;derval)
  }

i.dataFormat:{[x0]
  $[99h=type x0;value x0;
    0h >type x0;enlist x0;
    x0
  ]
  }
    

// minimize output of function by optimizing input variables
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/optimize.py#L1058
/* func    = function to be minimized
/* d    = dictionary of variables and arguments for func
/*        (keys are `xk and `arg-`arg only included if f takes other arguments than xk) 
/. r    > returns dictionary with optimal variables, corresponding function value and 
/.        number of iterations it took
newoptimize:{[func;x0;args;params]
  // update the default behaviour of the parameters
  params:updDefault[params];
  // format x0 based on input type
  x0:i.dataFormat[x0];
  // Evaluate the function at the starting point
  f0:funcEval[func;x0;args];
  // Calculate the starting gradient
  gk:i.grad[func;x0;args;f0;params`geps];
  // Initialize Hessian matrix as identity matrix
  hess:.ml.eye count x0;
  // set initial step guess i.e. the step before f0
  prev_fk:f0+sqrt[sum gk*gk]%2;
  gradNorm:vecNorm[gk;params`norm];
  optimKeys:`xk`fk`prev_fk`gk`prev_xk`hess`gnorm`I`idx;
  optimVals:(x0;f0;prev_fk;gk;0n;hess;gradNorm;hess;0);
  optimDict:optimKeys!optimVals;
  // Run optimization until one of the stopping conditions is met
  optimDict:stopCondition[;params]i.optimFunction[func;;args;params]/optimDict;
  returnKeys:`xVals`funcRet`numIter;
  // if function returned due to a null xVal or the new value being worse than the previous
  //  value then return the k-1 value
  returnVals:$[(optimDict[`fk]<optimDict`prev_fk) & (not any null optimDict`xk);
    optimDict`xk`fk`idx;
    optimDict`prev_xk`prev_fk`idx
    ];
  returnKeys!returnVals
  }

//optimimal func until gradient tolerance is reached
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/optimize.py#L1131
/* f        = function to be minimized
/* d        = dictionary of variables to be updated after each iterarion
/* fnc_keys = keys of variables passed to function 
/. r        > returns dictionary of optim variables and gradients at the end of each iteration 
i.optimFunction:{[f;d;args;params] 
  // calculate search direction
  pk:neg mmu[d`hess;d`gk];
  // line search func to be inserted to get alpha
  wolfe:i.wolfeSearch[;;;pk;f;;args;params]. d`fk`prev_fk`gk`xk;
  //old f_val goes to previous val
  d[`prev_fk]:d`fk;
  // update values from line search
  alpha   :wolfe 0;
  d[`fk]:wolfe 1;
  gnew    :wolfe 2;
  // define prev x_val
  d[`prev_xk]:d`xk;
  // update x values
  d[`xk]:d[`prev_xk]+alpha*pk;
  sk:d[`xk]-d`prev_xk;
  // if null gnew, then get gradient of new x value
  if[any null gnew;gnew:i.grad[f;d`xk;args;d`fk;params`geps]];
  // subtract new gradients
  yk:gnew-d`gk;
  d[`gk]:gnew;
  // get new norm of gradient
  d[`gnorm]:sqrt sum abs[d`gk]xexp 2;
  // calculate new hessian matrix for next iteration 
  rhok:1%mmu[yk;sk];
  A1:d[`I] - sk*\:yk*rhok;
  A2:d[`I] - yk*\:sk*rhok;
  d[`hess]:mmu[A1;mmu[d[`hess];A2]]+rhok*(sk*/:sk);
  // if returning infinite values end loop
  if[0w in abs d`xk;d[`gnorm`fk]:(0n;0w)];
  d[`idx]+:1;
  if[params`display;show d;-1"";];
  d
  }

// wolfeSearch func 
//   naming convention for dictionary keys is taken in this instance from the
//   python implementation outlined at 
//   https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L193
/* fk      = return value of f
/* prev_fk = prev value of f value before old
/* gk        = gradient values
/* pk        = search direction
/* d         = dict of x variables and arguments
/. r         > returns dictionary of new alpha, fk and derivative
i.wolfeSearch:{[fk;prev_fk;gk;pk;func;xk;args;params]
  phiFunc   :i.phi[func;pk;;xk;args];
  derphiFunc:i.derphi[func;params`geps;pk;;xk;args;fk];
  // initial Wolfe conditions
  wolfeDict:`idx`alpha0`phi0`phi_a0!(0;0;fk;fk);
  // calculate the derivative at that phi0
  derphi0:gk mmu pk;
  wolfeDict[`derphi_a0`derphi0]:2#derphi0;
  // the new alpha value should be between 0 and 1
  alphaval:1.01*2*(fk - prev_fk)%derphi0;
  // new alpha
  wolfeDict[`alpha1]:$[(alphaval>0) & (alphaval<1);alphaval;1];
  // phi value at alpha1
  wolfeDict[`phi_a1]:phiFunc wolfeDict`alpha1;
  // set up phi0 and derphi0
  // repeat until wolfe criteria is reached or max iteration
  // to get new alpha, phi and derphi values
  upd_d:{x[`idx]<y}[;10]i.scalarWolfe[derphiFunc;phiFunc;pk;params]/wolfeDict;
  // if the line search did not converge, use last alpha , phi and derphi
  $[not any null raze upd_d[`alpha_star`phi_star`derphi_star];
    upd_d[`alpha_star`phi_star`derphi_star];
    upd_d[`alpha1`phi_a1`derphi_a0_fin]]
 }

// Scalar search function search for Wolfe conditions
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L338
// This functions defines what are the "brackets" in between which the step function can be found.
// when optimal "bracket" is found, zoom in on area to find optimal
/* derphi_fnc = derivative fnc
/* phi_fnc    = phi function 
/* d          = dictionary with Wolfe values
/. r          > returns dictionary of new alpha, fk and derivative
i.scalarWolfe:{[derphiFunc;phiFunc;pk;params;d]
  // set up zoom function constant params
  zoom_setup:i.zoom_fnc[derphiFunc;phiFunc;;;params]. d`phi0`derphi0;
  // if criteria 1, zoom and break loop
  if[i.wolfeCriteria1[d;params];
    d[`idx]:0w;
    d[i.new_zoom]:zoom_setup d`alpha0`alpha1`phi_a0`phi_a1`derphi_a0;
    :d
    ];
  // calculate the derivative
  derphiCalc:derphiFunc d`alpha1;
  // update the new derivative fnc
  d[`derphi_a1]:derphiCalc`derval;
  // if criteria 2, then use current values for star values and break loop 
  $[i.wolfeCriteria2[d;params];
    [d[`alpha_star]:d`alpha1;
     d[`phi_star]:d`phi_a1;
     d[`derphi_star]:derphiCalc`grad;
     d[`idx]:0w;
     :d];
    // if criteria 3, zoom and stop loop
    0<=d`derphi_a1;
    [d[`idx]:0w;
     d[i.new_zoom]:zoom_setup d[`alpha1`alpha0`phi_a1`phi_a0`derphi_a1]];
    // update dictionary and repeat process until criteria is met
    [d[`alpha0]:d`alpha1;
     // new alpha has to be increased so just double
     d[`alpha1]:2*d`alpha1;
     d[`phi_a0]:d`phi_a1;
     d[`phi_a1]:phi_fnc[d`alpha1];
     d[`derphi_a0]:d`derphi_a1;
     d[`derphi_a0_fin]:derphiCalc`grad;
     d[`idx]+:1]
  ];
  d
  }

// Zoom in on "bracketed" area and find optimal step and fk for next step
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L537
/* derphi_fnc = derivative fnc
/* phi_fnc    = phi function 
/* phi0       = old value of f (fk)
/* derphi0    = inital derivative (grad * step direction)
/* lst        = list of hi and low calues for phi and deriv
/. r          > returns new alpha, fk and derivative
i.zoom_fnc:{[derphi_fnc;phi_fnc;phi0;derphi0;params;lst]
  d:i.z_dict!lst,phi0;
  d[`idx`a_rec]:2#0f;
  zoom_d:{x[`idx]<y}[;10]i.zoom[derphi_fnc;phi_fnc;phi0;derphi0;params]/d;
  // if zoom did not converge, set to null
  $[count star:zoom_d[i.new_zoom];star;3#0N]
  }

/^same as above, run below until criteria is met or max iterations is reached
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L556
i.zoom:{[derphi_fnc;phi_fnc;phi0;derphi0;params;d]
  // define high and low values
  dalpha:d[`a_hi]-d[`a_lo];
  h_l:`high`low!$[dalpha>0;d[`a_hi`a_lo];d[`a_lo`a_hi]];
  // Get cubic min
  fnd_min:i.cubicmin . d`a_lo`phi_lo`derphi_lo`a_hi`phi_hi`a_rec`phi_rec;
  // cubic interpolant chk
  cchk:dalpha*0.2;
  // if the result is too close to the end point points then use quadratic min 
  if[i.quadCriteria[fnd_min;h_l;cchk];
    // quadratic interpolant chk
    qchk:0.1*dalpha;
    fnd_min:i.quadmin . d`a_lo`phi_lo`derphi_lo`a_hi`phi_hi
  ];
  // update new values depending on fnd_min
  phi_min:phi_fnc[fnd_min];
  //first condition, update and continue loop
  if[i.zoomCriteria1[phi0;derphi0;phi_min;fnd_min;d;params];
    [d[`idx]+:1;
    d[i.upd_zoom1]:d[`phi_hi`a_hi],fnd_min,phi_min;:d]
  ];
  // calculate the derivative of the min
  derphi_min:derphi_fnc[fnd_min];
  // second scenario, create new features and end the loop
  $[i.zoomCriteria2[derphi0;derphi_min;params];
    [d[`idx]:0w;
     d:d,i.new_zoom!fnd_min,phi_min,enlist derphi_min`grad];
    // third condition, then update and continue in loop
    i.zoomCriteria3[derphi_min;dalpha];
    [d[`idx]+:1;
     d[i.upd_zoom1,i.upd_zoom2]:d[`phi_hi`a_hi`a_lo`phi_lo],fnd_min,phi_min,derphi_min`derval];
    //if no condition was satisfied update values and continue on loop
    [d[`idx]+:1;
     d[i.upd_zoom3,i.upd_zoom2]:d[`phi_lo`a_lo],fnd_min,phi_min,derphi_min[`derval]]
  ];
  // return updated dictionary
  d
  }


// Potentially can use lsq for this too? https://code.kx.com/q/ref/lsq/#polynomial-fitting

// minimizer of a cubic polynomial 
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L482
// minimize cubic function that goes through (a,fa), (b,fb), (c,fc)
// with derivative at a of fpa
// f(x) = A *(x-a)^3 + B*(x-a)^2 + C*(x-a) + D
// C = fpa
i.cubicmin:{[a;fa;fpa;b;fb;c;fc]
  db:b-a;
  dc:c-a;
  denom:(db*dc)xexp 2*(db-dc);
  d1:2 2#0f;
  d1[0]:(1 -1)*xexp[;2]each(db;dc);
  d1[1]:(-1 1)*xexp[;3]each(dc;db);
  AB:d1 mmu(fb-fa-fpa*db;fc-fa-fpa*dc);
  AB%:denom;
  radical:AB[1]*AB[1]-3*AB[0]*fpa;
  a+(neg[AB[1]]+sqrt(radical))%(3*AB[0])}

// get quadratic min
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L516
// get minimum of quad function that goes through points (a,fa), (b,fb)
// with derivative at point a of fpa
// f(x) = A*(x-a)^2 + B*(x-a) + C
// B = fpa
// C = fa
i.quadmin:{[a;fa;fpa;b;fb]
  db:b-a;
  B:(fb-fa-fpa*db)%(db*db);
  a-fpa%(2*B)
  }

// gradient func
/* fn   = function to be minimized
/* xk   = guess of value x at iteration k 
/* dict = dict of input values `xk`args!(x;args)
/* eps  = epsilon value
/* fk   = initial function evaluation completed 
/. r   > returns the gradient of each x variable
i.grad:{[func;xk;args;fk;eps]
  i.gradEval[fk;func;xk;args;eps]each til count xk
  }

// Criteria functionsi https://www.csie.ntu.edu.tw/~r97002/temp/num_optimization.pdf pg 60
i.wolfeCriteria1:{[d;params]
  (d[`phi_a1]>d[`phi0]+params[`c1]*d[`alpha1]*d[`derphi0])|((d[`phi_a1]>=d[`phi_a0])&(d[`idx]>1))
  }
i.wolfeCriteria2:{[d;params]
  neg[params[`c2]*d[`derphi0]]>=abs d`derphi_a1
  }

// Zoom criteria function https://www.csie.ntu.edu.tw/~r97002/temp/num_optimization.pdf pg 61
i.quadCriteria :{[fnd_min;h_l;cchk]
  (fnd_min>h_l[`low]-cchk)|fnd_min<h_l[`high]+cchk
  }
i.zoomCriteria1:{[phi0;derphi0;phi_min;fnd_min;d;params]
  (phi_min>phi0+fnd_min*derphi0*params`c1)|phi_min>=d[`phi_lo]
  }
i.zoomCriteria2:{[derphi0;derphi_min;params]
  abs[derphi_min`derval]<=neg derphi0*params`c2
  }
i.zoomCriteria3:{[derphi_min;dalpha]
  0 <= derphi_min[`derval]*dalpha
  }
   
   
// Zoom dictionary 

//input keys of zoom dictionary
i.z_dict:`a_lo`a_hi`phi_lo`phi_hi`derphi_lo`phi_rec;
// keys to be updated in zoom each iteration
i.upd_zoom1:`phi_rec`a_rec`a_hi`phi_hi;
// extra keys that have to be updated in some scenarios
i.upd_zoom2:`a_lo`phi_lo`derphi_lo;
i.upd_zoom3:`phi_rec`a_rec
// final updated keys to be used
i.new_zoom:`alpha_star`phi_star`derphi_star;
