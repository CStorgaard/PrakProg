# Berrut B₁ Rational Interpolation  
**Christian Storgaard’s Exam Project**  
Practical Programming with Numerical Methods

---

## 1. Introduction  
This project implements Berrut’s barycentric rational interpolants algorithm.  
- **First form** (\(B_1\)): reproduces constants exactly.  
- **Second form** (\(B_2\)): reproduces linear polynomials exactly.

## 2. Implementation Details  
- **Source**: `Berrut.cs`  
- **Classes**:  
  - `BerrutInterpolator(double[] nodes, double[] values, bool linearReproduction)`  
  - `Evaluate(double x)` returns the interpolated value.  
- **Weights**:  
  - First form:  
    \[
      w_j = (-1)^j, j = 0,...,n
    \]
  - Second form:  
      w_0 =  1/2
      w_n = (-1)^n{2} 
      w_j = (-1)^j, 1≤j≤n-1.

## 3. Files and running
Berrut.cs - Implementation of the interpolator, calculates weights, and evaluates values.
main.cs - sets up nodes & values for linear + sine function, perfoms interpolation, writes two data blocks out to Out.txt.
main_runge.cs - sets up nodes & values for the runge function, writes two data blocks out to Out_runge.txt.
main_time.cs - sets up different number of nodes for the linear + sine function and writes computation times to Out_times.txt.
main_error.cs - Calculated the absolute error between the two interpolators and the true function f(x) = x + 0.2sin(5x).
interpolation.gp - plots block 0 (nodes) and block 1 (r₁ & r₂).
runge.gp - plots block 0 (nodes) and block 1 (r₁ & r₂).
timing.gp - plots time in ns as a function of number of nodes to check computation time of O(n).
error.gp - plots the error between the interpolated values and real function f(x) = x + 0.2sin(5x) values.

Use make all, to make all relevant files.

## 4. Results
See interpolation.svg, runge.svg, time.svg and error.svg

Interpolation.svg:
Data (purple dots) is f(x)= x + 0.2 * sin(5*x) sampled at 9 equispaced nodes in [0,2π].
Red (r₁) is first form (constant exact), shows slight slope/offset bias.
Blue (r₂): second form (linear‐exact), sits exactly on the y=x baseline and adds the sine‐wiggles.

runge.svg
Data (purple dots) is the runge function f(x) = 1 / (1 + 25 * x**2) is sampled at 11 equidistant nodes from [-1,1].
For both interpolants, there are no significant oscillations at ±1, unlike higher degree polynomial fits, runge phenomenon is not present.

time.svg
Shows the O(n) relation between number of nodes placed and the evaluation time in ns.

error.svg
Shows the error between the function f(x) = x + 0.2sin(5x) and the interpolated values, gaps can be seen at node values as we are in log space and error there is 0.

## 5. Complexity and Properties
Precompute weights: O(n)

Evaluation cost: O(n) per point

Exactness: constants (first), lines (second).