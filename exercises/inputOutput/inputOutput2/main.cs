using static System.Console;
using static System.Math;
using System;


public class Program {
public static void Main(string[] args) {
char[] split_delimiters = {' ','\t','\n'};
var split_options = StringSplitOptions.RemoveEmptyEntries;
for( string line = ReadLine(); line != null; line = ReadLine() ){
	var numbers = line.Split(split_delimiters,split_options);
    WriteLine($"x \tsin(x) \tcos(x)");
	foreach(var number in numbers){
		double x = double.Parse(number);
		Error.WriteLine($"{x:F1} \t{Sin(x):F3} \t{Cos(x):F3}");
                }
        }
	}
}