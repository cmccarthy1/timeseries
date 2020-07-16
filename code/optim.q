// optimization function based off scipy https://docs.scipy.org/doc/scipy/reference/optimize.minimize-bfgs.html#optimize-minimize-bfgs
// Derivation of Formula http://www.bioinfo.org.cn/~wangchao/maa/Numerical_Optimization.pdf

// load in libraries
\l ml/ml.q
.ml.loadfile`:init.q

// Parameter for Armijo condition rule 
// and curvature condition rule respectively
c1:1e-4;c2:0.9;

// epsilon value for gradient fnc
// hardcoded for the moment, 
// might make it changable by the user in the future
eps:1.49e-8;

// minimize output of function by optimizing input variables
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/optimize.py#L1058
/* f    = function to be minimized
/* d    = dictionary of variables and arguments for func(keys are `xk and `arg-`arg only included if f takes other arguments than xk) 
/. r    > returns dictionary with optimimal variables and corresponding function value and how many iterations it took
optimize:{[f;d]
 // start value using start params
 fval:f . value d;
 // starting gradient
 gk:i.grad[f;d;eps];
 // Inital Hessan matrix is identity matrix
 hess:.ml.eye count d`xk;
 // set initial step guess (step before fval)
 prev_fval:fval+(sqrt sum gk*gk)%2;
 //norm of derivative
 gnorm:sqrt sum abs[gk]xexp 2;
 // keys and values to be passed to optimal function
 optim_keys:`fval`prev_fval`gk`prev_xk`hess`gnorm`I`i;
 optim_vals:(fval;prev_fval;gk;0n;hess;gnorm;hess;0);
 optim_dict:d,optim_keys!optim_vals;
 // run optim_fnc until grad tolerance is reached 
 optim_d:{x[`gnorm]>y}[;1e-4]i.optim_fnc[f;;key d]/optim_dict;
 // return the optimised parameters, function return value, and number of iterations it took
 opt_keys:`xk`fval`nit;
 // if increased fval, or xk is null then use previous value
 opt_vals:$[(optim_d[`fval]<optim_d[`prev_fval])&not any null optim_d[`xk];
           optim_d[`xk`fval`i];optim_d[`prev_xk`prev_fval`i]];
 opt_keys!opt_vals}

//optimimal func until gradient tolerance is reached
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/optimize.py#L1131
/* f        = function to be minimized
/* d        = dictionary of variables to be updated each iterarion
/* fnc_keys = keys of variables passed to function 
/. r        > returns dictionary of optim variables and gradients at the end of each iteration 
i.optim_fnc:{[f;d;fnc_keys] 
 // search direction
 pk:neg mmu[;]. d[`hess`gk];
 // line search func to be inserted to get alpha
 wolfe:i.wolfe_search[d`fval;d`prev_fval;d`gk;pk;f;fnc_keys!d[fnc_keys]];
 //old f_val goes to previous val
 d[`prev_fval]:d`fval;
 // update values from line search
 alpha:wolfe[0];
 d[`fval]:wolfe[1];
 gnew:wolfe[2];
 // define prev x_val
 d[`prev_xk]:d[`xk];
 // update x values
 d[`xk]:d[`prev_xk]+alpha*pk;
 sk:d[`xk]-d`prev_xk;
 // if null gnew, then get gradient of new x value
 if[any null gnew;gnew:i.grad[f;fnc_keys!d[fnc_keys];eps]];
 // subtract new gradients
 yk:gnew-d`gk;
 d[`gk]:gnew;
 // get new norm of gradient
 d[`gnorm]:sqrt sum abs[d[`gk]]xexp 2;
 // calculate new hessian matrix for next iteration 
 rhok:1%mmu[yk;sk];
 A1:d[`I] - sk*\:yk*rhok;
 A2:d[`I] - yk*\:sk*rhok;
 d[`hess]:mmu[A1;mmu[d[`hess];A2]]+rhok*(sk*/:sk);
 // if returning infinitiy values end loop
 if[0w in abs[d`xk];d[`gnorm`fval]:(0n;0w)];
 d[`i]+:1;d};

// wolfe_search func 
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L193
/* fval      = return value of f
/* prev_fval = prev value of f value before old
/* gk        = gradient values
/* pk        = search direction
/* d         = dict of x variables and arguments
/. r         > returns dictionary of new alpha, fval and derivative
i.wolfe_search:{[fval;prev_fval;gk;pk;f;d]
 // set up phi and derphi fncs
 phi_fnc:i.phi[f;pk;;d];
 derphi_fnc:i.derphi[f;eps;pk;;d];
 // set up for wolfe func
 wolfe_d:`alpha0`alpha1`phi_a0`phi_a1`derphi_a0`i!();
 // prev alpha initally set to 0
 wolfe_d[`alpha0]:0;
 // get the derivative at that point
 derphi0:gk mmu pk;
 // the new alpha value should be between 0 and 1
 alphaval:1.01*2*(fval - prev_fval)%derphi0;
 // new alpha
 wolfe_d[`alpha1]:$[(alphaval>0)&(alphaval<1);alphaval;1];
 // phi value at alpha1
 wolfe_d[`phi_a1]:phi_fnc[wolfe_d[`alpha1]];
 wolfe_d[`phi_a0]:fval;
 // derphi initial value
 wolfe_d[`derphi_a0]:derphi0;
 wolfe_d[`i]:1;
 // set up phi0 and derphi0
 wolfe_d[`phi0`derphi0]:fval,derphi0;
 // repeat until wolfe criteria is reached or max iteration
 // to get new alpha, phi and derphi values
 upd_d:{x[`i]<y}[;10]i.scalar_wolfe[derphi_fnc;phi_fnc;pk]/wolfe_d;
 // if the line search did not converge, use last alpha , phi and derphi
 $[not any null raze upd_d[`alpha_star`phi_star`derphi_star];
    upd_d[`alpha_star`phi_star`derphi_star];upd_d[`alpha1`phi_a1`derphi_a0_fin]]
 }

// Scalar search function search for Wolfe conditions
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L338
// This functions defines what are the "brackets" in between which the step function can be found.
// when optimal "bracket" is found, zoom in on area to find optimal
/* derphi_fnc = derivative fnc
/* phi_fnc    = phi function 
/* d          = dictionary with Wolfe values
/. r          > returns dictionary of new alpha, fval and derivative
i.scalar_wolfe:{[derphi_fnc;phi_fnc;pk;d]
  // set up zoom function constant params
  zoom_setup:i.zoom_fnc[derphi_fnc;phi_fnc;;]. d`phi0`derphi0;
  // if criteria 1, zoom and break loop
  if[i.wolfe_crit_1[d];
   [d[`i]:0w;
    d[i.new_zoom]:zoom_setup d[`alpha0`alpha1`phi_a0`phi_a1`derphi_a0];:d]];
    // calculate the derivative
    derphi_calc:derphi_fnc d`alpha1;
   // update the new derivative fnc
    d[`derphi_a1]:derphi_calc`derval;
    // if criteria 2, then use current values for star values and break loop 
     $[i.wolfe_crit_2[d];
     [d[`alpha_star]:d[`alpha1];
      d[`phi_star]:d[`phi_a1];
      d[`derphi_star]:derphi_calc`grad;
      d[`i]:0w;d];
   // if criteria 3, zoom and stop loop
     0<=d[`derphi_a1];
    [d[`i]:0w;
     d[i.new_zoom]:zoom_setup d[`alpha1`alpha0`phi_a1`phi_a0`derphi_a1]];
   // update dictionary and repeat process until criteria is met
   [d[`alpha0]:d[`alpha1];
    // new alpha has to be increased so just double
    d[`alpha1]:2*d[`alpha1];
    d[`phi_a0]:d[`phi_a1];
    d[`phi_a1]:phi_fnc[d`alpha1];
    d[`derphi_a0]:d`derphi_a1;
    d[`derphi_a0_fin]:derphi_calc`grad;
    d[`i]+:1]];d
 }

// Zoom in on "bracketed" area and find optimal step and fval for next step
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L537
/* derphi_fnc = derivative fnc
/* phi_fnc    = phi function 
/* phi0       = old value of f (fval)
/* derphi0    = inital derivative (grad * step direction)
/* lst        = list of hi and low calues for phi and deriv
/. r          > returns new alpha, fval and derivative
i.zoom_fnc:{[derphi_fnc;phi_fnc;phi0;derphi0;lst]
  d:i.z_dict!lst,phi0;
  d[`i`a_rec]:2#0f;
  zoom_d:{x[`i]<y}[;10]i.zoom[derphi_fnc;phi_fnc;phi0;derphi0]/d;
  // if zoom did not converge, set to null
 $[count star:zoom_d[i.new_zoom];star;3#0N]
  }

/^same as above, run below until criteria is met or max iterations is reached
// https://github.com/scipy/scipy/blob/v1.5.0/scipy/optimize/linesearch.py#L556
i.zoom:{[derphi_fnc;phi_fnc;phi0;derphi0;d]
 // define high and low values
 dalpha:d[`a_hi]-d[`a_lo];
 h_l:`high`low!$[dalpha>0;d[`a_hi`a_lo];d[`a_lo`a_hi]];
 // Get cubic min
 fnd_min:i.cubicmin . d`a_lo`phi_lo`derphi_lo`a_hi`phi_hi`a_rec`phi_rec;
 // cubic interpolant chk
 cchk:dalpha*0.2;
 // if the result is too close to the end point points then use quadratic min 
 if[i.quad_crit[fnd_min;h_l;cchk];
   // quadratic interpolant chk
   qchk:0.1*dalpha;
   fnd_min:i.quadmin . d`a_lo`phi_lo`derphi_lo`a_hi`phi_hi];
  // update new values depending on fnd_min
  phi_min:phi_fnc[fnd_min];
  //first condition, update and continue loop
  if[i.zoom_crit_1[phi0;derphi0;phi_min;fnd_min;d];
    [d[`i]+:1;
     d[i.upd_zoom1]:d[`phi_hi`a_hi],fnd_min,phi_min;:d]];
  // calculate the derivative of the min
   derphi_min:derphi_fnc[fnd_min];
  // second scenari, create new features and end the loop
  $[i.zoom_crit_2[derphi0;derphi_min];
   [d[`i]:0w;
   d:d,i.new_zoom!fnd_min,phi_min,enlist[d_grad:derphi_min[`grad]]];
  // third condition, then update and continue in loop
  i.zoom_crit_3[derphi_min;dalpha];
    [d[`i]+:1;
    d[i.upd_zoom1,i.upd_zoom2]:d[`phi_hi`a_hi`a_lo`phi_lo],fnd_min,phi_min,derphi_min`derval];
  //if no condition was satisfied update values and continue on loop
  [d[`i]+:1;
  d[i.upd_zoom3,i.upd_zoom2]:d[`phi_lo`a_lo],fnd_min,phi_min,derphi_min[`derval]]];
  // return updated dictionary
  d
 }

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
 d1[0;0]:dc xexp 2;
 d1[0;1]:neg db xexp 2;
 d1[1;0]:neg dc xexp 3;
 d1[1;1]:db xexp 3;
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
/* f   = function to be minimized
/* inp = dict of input values (x;args)
/* eps = epsilon value
/. r   > returns the gradient of each x variable
i.grad:{[f;inp;eps]
   {[f;inp;eps;i]
   // inital f value
   f0:f . inp;
   // increase x@i by eps
   inp[0;i]:inp[0;i]+eps;
   // increase in new f[x] divided by increase in x value
   ((f . inp)-f0)%eps
   }[f;value inp;eps]each til count inp`xk}

// function with step value applied
/* f     = function to be minimized
/* pk    = step direction
/* alpha = step applied
/* d     > dictionary of x variables and argument
/. r - returns value of function with step value applied
i.phi:{[f;pk;alpha;d]
  d[`xk]+:alpha*pk;
  f . value d}

// derivative of function with step value included
/* f     = function to be minimized
/* pk    = step direction
/* alpha = step applied
/* d     = dictionary of x variables and argument
/. r     > returns dictionary of gradiant and scalar derivative
i.derphi:{[f;eps;pk;alpha;d]
   // increase xk 
   d[`xk]+:alpha*pk;
   // get gradient
   gval:i.grad[f;d;eps];
   derval:gval mmu pk;
   `grad`derval!(gval;derval)}

// Criteria functionsi https://www.csie.ntu.edu.tw/~r97002/temp/num_optimization.pdf pg 60
i.wolfe_crit_1:{[d](d[`phi_a1]>d[`phi0]+c1*d[`alpha1]*d[`derphi0])|((d[`phi_a1]>=d[`phi_a0])&(d[`i]>1))};
i.wolfe_crit_2:{[d]neg[c2*d[`derphi0]]>=abs d`derphi_a1}

// Zoom criteria function https://www.csie.ntu.edu.tw/~r97002/temp/num_optimization.pdf pg 61
i.quad_crit:  {[fnd_min;h_l;cchk](fnd_min>h_l[`low]-cchk)|fnd_min<h_l[`high]+cchk}
i.zoom_crit_1:{[phi0;derphi0;phi_min;fnd_min;d](phi_min>phi0+c1*fnd_min*derphi0)|phi_min>=d[`phi_lo]};
i.zoom_crit_2:{[derphi0;derphi_min]abs(derphi_min`derval)<=neg[c2*derphi0]}
i.zoom_crit_3:{[derphi_min;dalpha](derphi_min[`derval]*dalpha)>=0}
   
   
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
