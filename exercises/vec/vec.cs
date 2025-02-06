using System;


/// <summary>
/// Represents a 3D vector and provides various vector operations.
/// </summary>
public class Vec
{
    public double x, y, z; // Components of the vector

    // Default constructor initializes vector to (0,0,0)
    public Vec() { x = y = z = 0; }

    // Parameterized constructor
    public Vec(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }

    // Scalar multiplication
    public static Vec operator *(Vec v, double c) => new Vec(c * v.x, c * v.y, c * v.z);
    public static Vec operator *(double c, Vec v) => v * c;

    // Scalar division
    public static Vec operator /(Vec v, double c) => new Vec(v.x / c, v.y / c, v.z / c);

    // Vector addition
    public static Vec operator +(Vec u, Vec v) => new Vec(u.x + v.x, u.y + v.y, u.z + v.z);

    // Vector negation
    public static Vec operator -(Vec u) => new Vec(-u.x, -u.y, -u.z);

    // Vector subtraction
    public static Vec operator -(Vec u, Vec v) => new Vec(u.x - v.x, u.y - v.y, u.z - v.z);

    /// <summary>
    /// Computes the dot product of this vector with another vector.
    /// </summary>
    public double Dot(Vec other) => this.x * other.x + this.y * other.y + this.z * other.z;

    /// <summary>
    /// Computes the dot product of two vectors.
    /// </summary>
    public static double Dot(Vec v, Vec w) => v.x * w.x + v.y * w.y + v.z * w.z;

    /// <summary>
    /// Approximate equality comparison for doubles.
    /// </summary>
    public static bool Approx(double a, double b, double acc = 1e-9, double eps = 1e-9)
    {
        if (Math.Abs(a - b) < acc) return true;
        if (Math.Abs(a - b) < (Math.Abs(a) + Math.Abs(b)) * eps) return true;
        return false;
    }

    /// <summary>
    /// Approximate equality comparison for vectors.
    /// </summary>
    public bool Approx(Vec other) => Approx(this.x, other.x) && Approx(this.y, other.y) && Approx(this.z, other.z);
    public static bool Approx(Vec u, Vec v) => u.Approx(v);

    /// <summary>
    /// Returns a string representation of the vector.
    /// </summary>
    public override string ToString() => $"{x} {y} {z}";

        /// <summary>
    /// Prints the vector to the console.
    /// </summary>
    public void Print(string s = "")
{
    Console.Write(s);
    Console.WriteLine($"{x} {y} {z}");
}
}


