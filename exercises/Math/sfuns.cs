using System;

public class sfuns{

// Stirling Approximation of Gamma-function

public static double fgamma(double x){
///single precision gamma function (formula from Wikipedia)
if(x<0)return Math.PI/Math.Sin(Math.PI*x)/fgamma(1-x); // Euler's reflection formula
if(x<9)return fgamma(x+1)/x; // Recurrence relation
double lnfgamma=x*Math.Log(x+1/(12*x-1/x/10))-x+Math.Log(2*Math.PI/x)/2;
return Math.Exp(lnfgamma);
}

public static double lnfgamma(double x){
///single precision gamma function (formula from Wikipedia)
if(x <= 0) return double.NaN;
if(x < 9) return lnfgamma(x+1) - Math.Log(x);
double res=x*Math.Log(x+1/(12*x-1/x/10))-x+Math.Log(2*Math.PI/x)/2;
return res;
}

}