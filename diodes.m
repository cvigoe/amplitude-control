function [R3, R4, r1, r2, v_d1, v_d2] = diodes(R3_0, R4_0, r1_0, r2_0, v_d1_0, v_d2_0)
%[t1, t2, t3, t4, t5, t6, t7, t8] = diodes(R3_0, R4_0, r1_0, r2_0, v_d1_0, v_d2_0): 
%function to approximate the vector of unknowns X for the non-linear
%amplitude control circuit implemented in the Alternative Phase Shift 
%Oscillator circuit of the Oscillator Design Lab for EEEN30120
%using the multi-dimensional Newton-Raphson method with
%initial estimates R3_0, R4_0, r1_0, r2_0, v_d1_0, v_d2_0
%
%Input R3_0, R4_0 = initial estimates for resistances of resistors (ohms)
%Input r1_0, r2_0 = initial estimates effective resistances of diodes (ohms)
%Input v_d1_0, v_d2_0 = initial estimates for diode voltages (V)
%Output R3, R4 = approximate resistances for R3 and R4 (ohms)
%Output r1, r2 = approximate effective resistances for diodes at bias point (ohms)
%Output v_d1, v_d2 = approximate diode voltages (V)

%   Version 1: created 09/03/2017. Author: Conor Igoe
%   This MATLAB function M-file is not flexible. It works for the
%   non-linear amplitude control circuit used in the Alternative Phase
%   Shift Oscillator circuit only. Note that as precision is not required,
%   scaling is not used in Newton-Raphson, limiting its maximum accuracy.
%
%   The limit on the number of allowable iterations and the acceptable 
%   tolerance for the Newton-Raphson (NR) algorithm are internally 
%   generated.

%% Check input and output arguments
if (nargin ~= 6), error('Incorrect number of input arguments.'); end
if (nargout ~= 6), error('Incorrect number of output arguments.'); end

%% Internal parameter ITERATION_LIMIT = maximum number of steps permitted.
% Internal parameter TOLERANCE = minimum acceptable value for |f_dash(x)|
%                                evaluated at the approximate roots.
% Internal parameter I_s = diode saturation current for use in Shockley
% Diode Equation
% Internal parameter V_T = thermal voltage for use in Diode Equation
% Internal parameter V_out = peak of required sinusoidal output

ITERATION_LIMIT = 10;
TOLERANCE = 10^-1;

I_s = 10e-12;
V_T = 26e-3;
V_out = 3.5;

%% System of equations to be used in NR
f1 = @(R3, R4, r1, r2, v_d1, v_d2) ( (I_s/V_T)*exp(v_d1/V_T) - 1/(r1));
f2 = @(R3, R4, r1, r2, v_d1, v_d2) ( (I_s/V_T)*exp(v_d2/V_T) - 1/(r2));
f3 = @(R3, R4, r1, r2, v_d1, v_d2) ( V_out*R4*r1*r2/(R3*(R4*(r1 + r2) + r1*r2) + R4*r1*r2)  - v_d1);
f4 = @(R3, R4, r1, r2, v_d1, v_d2) ( -V_out*R4*r1*r2/(R3*(R4*(r1 + r2) + r1*r2) + R4*r1*r2) - v_d2);
f5 = @(R3, R4, r1, r2, v_d1, v_d2) ( R3 + R4 - 68000 );
f6 = @(R3, R4, r1, r2, v_d1, v_d2) ( R3 + (R4*r1*r2/(R4*(r1 + r2) + r1*r2)) - 47500 );

%% Jacobian of non-linear system of equations to be used in NR
d_f1_R3 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f1_R4 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f1_r1 = @(R3, R4, r1, r2, v_d1, v_d2) ( 1/(r1*r1) );
d_f1_r2 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f1_v_d1 = @(R3, R4, r1, r2, v_d1, v_d2) ( I_s/(V_T*V_T)*exp(v_d1/V_T) );
d_f1_v_d2 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );

d_f2_R3 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f2_R4 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f2_r1 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f2_r2 = @(R3, R4, r1, r2, v_d1, v_d2) ( 1/(r2*r2) );
d_f2_v_d1 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f2_v_d2 = @(R3, R4, r1, r2, v_d1, v_d2) ( I_s/(V_T*V_T)*exp(v_d2/V_T) );

d_f3_R3 = @(R3, R4, r1, r2, v_d1, v_d2) ( -V_out*((R4*r1*r2/(R4*(r1 + r2) + r1*r2)))/((R3 + ((R4*r1*r2/(R4*(r1 + r2) + r1*r2))))*((R3 + ((R4*r1*r2/(R4*(r1 + r2) + r1*r2)))))) );
d_f3_R4 = @(R3, R4, r1, r2, v_d1, v_d2) ( V_out*(r1*r2/(r1 + r2))*(r1*r2/(r1 + r2))*R3 / (((r1*r2/(r1 + r2))*R3 + (r1*r2/(r1 + r2))*R4 + R3*R4)*((r1*r2/(r1 + r2))*R3 + (r1*r2/(r1 + r2))*R4 + R3*R4)) );
d_f3_r1 = @(R3, R4, r1, r2, v_d1, v_d2) ( V_out*(R4*r2/(R4 + r2))*(R4*r2/(R4 + r2))*R3 / (((R4*r2/(R4 + r2))*R3 + (R4*r2/(R4 + r2))*r1 + R3*r1)*((R4*r2/(R4 + r2))*R3 + (R4*r2/(R4 + r2))*r1 + R3*r1)) );
d_f3_r2 = @(R3, R4, r1, r2, v_d1, v_d2) ( V_out*(r1*R4/(r1 + R4))*(r1*R4/(r1 + R4))*R3 / (((r1*R4/(r1 + R4))*R3 + (r1*R4/(r1 + R4))*r2 + R3*r2)*((r1*R4/(r1 + R4))*R3 + (r1*R4/(r1 + R4))*r2 + R3*r2)) );
d_f3_v_d1 = @(R3, R4, r1, r2, v_d1, v_d2) ( -1 );
d_f3_v_d2 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );

d_f4_R3 = @(R3, R4, r1, r2, v_d1, v_d2) ( V_out*((R4*r1*r2/(R4*(r1 + r2) + r1*r2)))/((R3 + ((R4*r1*r2/(R4*(r1 + r2) + r1*r2))))*((R3 + ((R4*r1*r2/(R4*(r1 + r2) + r1*r2)))))) );
d_f4_R4 = @(R3, R4, r1, r2, v_d1, v_d2) ( -V_out*(r1*r2/(r1 + r2))*(r1*r2/(r1 + r2))*R3 / (((r1*r2/(r1 + r2))*R3 + (r1*r2/(r1 + r2))*R4 + R3*R4)*((r1*r2/(r1 + r2))*R3 + (r1*r2/(r1 + r2))*R4 + R3*R4)) );
d_f4_r1 = @(R3, R4, r1, r2, v_d1, v_d2) ( -V_out*(R4*r2/(R4 + r2))*(R4*r2/(R4 + r2))*R3 / (((R4*r2/(R4 + r2))*R3 + (R4*r2/(R4 + r2))*r1 + R3*r1)*((R4*r2/(R4 + r2))*R3 + (R4*r2/(R4 + r2))*r1 + R3*r1)) );
d_f4_r2 = @(R3, R4, r1, r2, v_d1, v_d2) ( -V_out*(r1*R4/(r1 + R4))*(r1*R4/(r1 + R4))*R3 / (((r1*R4/(r1 + R4))*R3 + (r1*R4/(r1 + R4))*r2 + R3*r2)*((r1*R4/(r1 + R4))*R3 + (r1*R4/(r1 + R4))*r2 + R3*r2)) );
d_f4_v_d1 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f4_v_d2 = @(R3, R4, r1, r2, v_d1, v_d2) ( -1 );

d_f5_R3 = @(R3, R4, r1, r2, v_d1, v_d2) ( 1 );
d_f5_R4 = @(R3, R4, r1, r2, v_d1, v_d2) ( 1 );
d_f5_r1 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f5_r2 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f5_v_d1 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f5_v_d2 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );

d_f6_R3 = @(R3, R4, r1, r2, v_d1, v_d2) ( 1 );
d_f6_R4 = @(R3, R4, r1, r2, v_d1, v_d2) ( ((r1*r2)/(r1 + r2))*((r1*r2)/(r1 + r2))/( (((r1*r2)/(r1 + r2)) + R4) * (((r1*r2)/(r1 + r2)) + R4)  ) );
d_f6_r1 = @(R3, R4, r1, r2, v_d1, v_d2) ( ((R4*r2)/(R4 + r2))*((R4*r2)/(R4 + r2))/( (((R4*r2)/(R4 + r2)) + r1) * (((R4*r2)/(R4 + r2)) + r1)  ) );
d_f6_r2 = @(R3, R4, r1, r2, v_d1, v_d2) ( ((r1*R4)/(r1 + R4))*((r1*R4)/(r1 + R4))/( (((r1*R4)/(r1 + R4)) + r2) * (((r1*R4)/(r1 + R4)) + r2)  ) );
d_f6_v_d1 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );
d_f6_v_d2 = @(R3, R4, r1, r2, v_d1, v_d2) ( 0 );

%% Newton-Raphson Algorithm
X = [R3_0; R4_0; r1_0; r2_0; v_d1_0; v_d2_0];       % Initialize initial guess
iterations = 0;                                     % Initialize iterations counter

for count = 1:ITERATION_LIMIT + 1
    if count == ITERATION_LIMIT + 1
        error('The Iteration limit has been reached and did not converge')
        break
    end
    
    % Update root estimate
    R3 = X(1);
    R4 = X(2);
    r1 = X(3);
    r2 = X(4);
    v_d1 = X(5);
    v_d2 = X(6);
    
    % Check if we have met the desired tolerance
    if abs(f1(R3, R4, r1, r2, v_d1, v_d2)) < TOLERANCE && abs(f2(R3, R4, r1, r2, v_d1, v_d2)) < TOLERANCE && abs(f3(R3, R4, r1, r2, v_d1, v_d2)) < TOLERANCE && abs(f4(R3, R4, r1, r2, v_d1, v_d2)) < TOLERANCE && abs(f5(R3, R4, r1, r2, v_d1, v_d2)) < TOLERANCE && abs(f6(R3, R4, r1, r2, v_d1, v_d2)) < TOLERANCE
        break
    end
    
    % Calculate f(x)
    f = [f1(R3, R4, r1, r2, v_d1, v_d2); f2(R3, R4, r1, r2, v_d1, v_d2); f3(R3, R4, r1, r2, v_d1, v_d2); f4(R3, R4, r1, r2, v_d1, v_d2); f5(R3, R4, r1, r2, v_d1, v_d2); f6(R3, R4, r1, r2, v_d1, v_d2)];
    
    % Calculate D_f_dash(x)
    D_f =      [[d_f1_R3(R3, R4, r1, r2, v_d1, v_d2) d_f1_R4(R3, R4, r1, r2, v_d1, v_d2) d_f1_r1(R3, R4, r1, r2, v_d1, v_d2) d_f1_r2(R3, R4, r1, r2, v_d1, v_d2) d_f1_v_d1(R3, R4, r1, r2, v_d1, v_d2) d_f1_v_d2(R3, R4, r1, r2, v_d1, v_d2)];
                [d_f2_R3(R3, R4, r1, r2, v_d1, v_d2) d_f2_R4(R3, R4, r1, r2, v_d1, v_d2) d_f2_r1(R3, R4, r1, r2, v_d1, v_d2) d_f2_r2(R3, R4, r1, r2, v_d1, v_d2) d_f2_v_d1(R3, R4, r1, r2, v_d1, v_d2) d_f2_v_d2(R3, R4, r1, r2, v_d1, v_d2)];
                [d_f3_R3(R3, R4, r1, r2, v_d1, v_d2) d_f3_R4(R3, R4, r1, r2, v_d1, v_d2) d_f3_r1(R3, R4, r1, r2, v_d1, v_d2) d_f3_r2(R3, R4, r1, r2, v_d1, v_d2) d_f3_v_d1(R3, R4, r1, r2, v_d1, v_d2) d_f3_v_d2(R3, R4, r1, r2, v_d1, v_d2)];
                [d_f4_R3(R3, R4, r1, r2, v_d1, v_d2) d_f4_R4(R3, R4, r1, r2, v_d1, v_d2) d_f4_r1(R3, R4, r1, r2, v_d1, v_d2) d_f4_r2(R3, R4, r1, r2, v_d1, v_d2) d_f4_v_d1(R3, R4, r1, r2, v_d1, v_d2) d_f4_v_d2(R3, R4, r1, r2, v_d1, v_d2)];
                [d_f5_R3(R3, R4, r1, r2, v_d1, v_d2) d_f5_R4(R3, R4, r1, r2, v_d1, v_d2) d_f5_r1(R3, R4, r1, r2, v_d1, v_d2) d_f5_r2(R3, R4, r1, r2, v_d1, v_d2) d_f5_v_d1(R3, R4, r1, r2, v_d1, v_d2) d_f5_v_d2(R3, R4, r1, r2, v_d1, v_d2)];
                [d_f6_R3(R3, R4, r1, r2, v_d1, v_d2) d_f6_R4(R3, R4, r1, r2, v_d1, v_d2) d_f6_r1(R3, R4, r1, r2, v_d1, v_d2) d_f6_r2(R3, R4, r1, r2, v_d1, v_d2) d_f6_v_d1(R3, R4, r1, r2, v_d1, v_d2) d_f6_v_d2(R3, R4, r1, r2, v_d1, v_d2)]];

    % Calculate next approximation
    X = X - (D_f)\f;
   
    iterations = iterations + 1;
end
end
