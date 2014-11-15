% Model: Basic RBC model with endogenous labor supply
%   
% Author: Dario Caldara, based on codes by Jesus Fernandez-Villaverde
%
% Purpose: Prepared for the course: "Tools for Nonlinear DSGE models 
% Last Update: 29.06.2010
%
%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

close all

%----------------------------------------------------------------
% 1. Endogenous variables
%----------------------------------------------------------------
var 

// Allocation variables 
k y c l i  

// Utility variables
u v

// Input prices
r w 

rf

// Shocks
z
;

%----------------------------------------------------------------
% 2. Exogenous variables
%----------------------------------------------------------------

varexo e;

%----------------------------------------------------------------
% 3. Parameters
%----------------------------------------------------------------

parameters 

// Utility function
nu beta gamma

// Technology 
alpha delta lambda sigma 


// Scaling factor
cte
;
%----------------------------------------------------------------
% 4. Calibration
%----------------------------------------------------------------

// Utility 
beta = 0.991;
nu =  0.3622;
gamma = 2; 
//gamma = 50; 

// Technology
alpha = 0.3;
delta = 0.0196;
lambda = 0.95;
// sigma = 0.007;
sigma = 0.021;

%----------------------------------------------------------------
% 5. Steady State
%----------------------------------------------------------------
A = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
B = (nu*(1-alpha)*A^(alpha))/((1-nu)*(A^(alpha-1)-delta));

l_ss = B/(A+B);
k_ss  = B*(1-l_ss);
c_ss = k_ss^alpha*l_ss^(1-alpha)-delta*k_ss;
w_ss = (1-alpha)*k_ss^alpha*l_ss^(-alpha);
r_ss  = 1+alpha*k_ss^(alpha-1)*l_ss^(1-alpha)-delta;

cte = (c_ss^nu*(1-l_ss)^(1-nu))^(-1);
u_ss = cte*c_ss^nu*(1-l_ss)^(1-nu);
v_ss = u_ss;

%----------------------------------------------------------------
% 6. Model
%----------------------------------------------------------------

//model(use_dll); 
model; 

    // 1. Utitlity
    u = cte*(c^nu*(1-l)^(1-nu));
    
    // 2. Production Function
    y = exp(z)*k(-1)^alpha*l^(1-alpha);
    
    // 3. Law of Motion for Productivity
    z = lambda*z(-1) + sigma*e;

    // 4. Law of Motion for Capital
    k = (1-delta)*k(-1) + i;
    
    // 5. Resource Constraint
    c + i = y;

    // 6. Euler Equation for Capital
    beta*((u(+1)/u)^(1-gamma))*(c/c(+1))*r(+1) = 1;

    // 7. Static Leisure-Consumption
    (1-nu)/nu*c/(1-l)= w;

    // 8. Return on Risky Asset
    r = 1+alpha*y/k(-1) - delta;

    // 9. Wage
    w = (1-alpha)*y/l;

    // 10. Value Function
    v = (1-beta)/(1-gamma)*u^(1-gamma) + beta*v(+1);

    // 11. Risk-free Rate
    beta*((u(+1)/u)^(1-gamma))*(c/c(+1)) = 1/rf;


end;

%----------------------------------------------------------------
% 7. Computation
%----------------------------------------------------------------

initval;
  l = l_ss;
  k = k_ss;
  c = c_ss;
  w = w_ss;
  r = r_ss;
  u = v_ss;
  v = v_ss;
  rf= 1/beta;

  z = 0;
  e = 0;

end;
    
shocks;
    var e = 1;
end;

steady;

set_dynare_seed(2)
stoch_simul(PRUNING, PERIODS = 100000, DROP = 1000, IRF = 0, ORDER = 2) c y i k rf r;