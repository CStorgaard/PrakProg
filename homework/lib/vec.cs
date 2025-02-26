using System;
using System.Text;

/// <summary>
/// Represents a vector of arbitrary length and provides various vector operations.
/// </summary>
public class Vec
{
    private double[] data;

    /// <summary>
    /// Gets the dimension (number of components) of the vector.
    /// </summary>
    public int Length => data.Length;

    /// <summary>
    /// Default constructor initializes a 3-dimensional zero vector.
    /// </summary>
    public Vec() : this(3) { }

    /// <summary>
    /// Constructs a vector with the specified dimension, initializing all components to zero.
    /// </summary>
    public Vec(int dimension)
    {
        if (dimension <= 0)
            throw new ArgumentException("Dimension must be positive.", nameof(dimension));
        data = new double[dimension];
    }

    /// <summary>
    /// Constructs a vector from an array of doubles.
    /// </summary>
    public Vec(double[] values)
    {
        if (values == null || values.Length == 0)
            throw new ArgumentException("Values cannot be null or empty.", nameof(values));
        data = new double[values.Length];
        Array.Copy(values, data, values.Length);
    }

    /// <summary>
    /// Constructs a 3-dimensional vector with the specified x, y, and z components.
    /// </summary>
    public Vec(double x, double y, double z) : this(new double[] { x, y, z }) { }

    /// <summary>
    /// Indexer to get or set the components of the vector.
    /// </summary>
    public double this[int index]
    {
        get 
        {
            if (index < 0 || index >= Length)
                throw new IndexOutOfRangeException();
            return data[index];
        }
        set 
        {
            if (index < 0 || index >= Length)
                throw new IndexOutOfRangeException();
            data[index] = value;
        }
    }

    /// <summary>
    /// Sets the element at the specified index to the given value.
    /// </summary>
    public void SetElement(int index, double value)
    {
        this[index] = value;
    }

    /// <summary>
    /// Computes the Euclidean norm (magnitude) of the vector.
    /// </summary>
    public double Norm()
    {
        double sum = 0;
        for (int i = 0; i < Length; i++)
            sum += data[i] * data[i];
        return Math.Sqrt(sum);
    }

    /// <summary>
    /// Multiplies a vector by a scalar.
    /// </summary>
    public static Vec operator *(Vec v, double c)
    {
        Vec result = new Vec(v.Length);
        for (int i = 0; i < v.Length; i++)
            result[i] = v[i] * c;
        return result;
    }

    /// <summary>
    /// Multiplies a scalar by a vector.
    /// </summary>
    public static Vec operator *(double c, Vec v) => v * c;

    /// <summary>
    /// Divides a vector by a scalar.
    /// </summary>
    public static Vec operator /(Vec v, double c)
    {
        if (c == 0)
            throw new DivideByZeroException("Division by zero is not allowed.");
        Vec result = new Vec(v.Length);
        for (int i = 0; i < v.Length; i++)
            result[i] = v[i] / c;
        return result;
    }

    /// <summary>
    /// Adds two vectors element-wise.
    /// </summary>
    public static Vec operator +(Vec u, Vec v)
    {
        if (u.Length != v.Length)
            throw new ArgumentException("Vectors must be of the same dimension for addition.");
        Vec result = new Vec(u.Length);
        for (int i = 0; i < u.Length; i++)
            result[i] = u[i] + v[i];
        return result;
    }

    /// <summary>
    /// Negates the vector.
    /// </summary>
    public static Vec operator -(Vec v) => v * -1;

    /// <summary>
    /// Subtracts one vector from another element-wise.
    /// </summary>
    public static Vec operator -(Vec u, Vec v)
    {
        if (u.Length != v.Length)
            throw new ArgumentException("Vectors must be of the same dimension for subtraction.");
        Vec result = new Vec(u.Length);
        for (int i = 0; i < u.Length; i++)
            result[i] = u[i] - v[i];
        return result;
    }

    /// <summary>
    /// Computes the dot product of this vector with another vector.
    /// </summary>
    public double Dot(Vec other)
    {
        if (other == null)
            throw new ArgumentNullException(nameof(other));
        if (this.Length != other.Length)
            throw new ArgumentException("Vectors must be of the same dimension for dot product.");
        double sum = 0;
        for (int i = 0; i < Length; i++)
            sum += this[i] * other[i];
        return sum;
    }

    /// <summary>
    /// Computes the dot product of two vectors.
    /// </summary>
    public static double Dot(Vec u, Vec v) => u.Dot(v);

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
    public bool Approx(Vec other)
    {
        if (other == null || this.Length != other.Length)
            return false;
        for (int i = 0; i < this.Length; i++)
        {
            if (!Approx(this[i], other[i]))
                return false;
        }
        return true;
    }

    public static bool Approx(Vec u, Vec v) => u.Approx(v);

    /// <summary>
    /// Returns a string representation of the vector.
    /// </summary>
    public override string ToString()
    {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < Length; i++)
        {
            sb.Append(this[i].ToString("G"));
            if (i < Length - 1)
                sb.Append(" ");
        }
        return sb.ToString();
    }

    /// <summary>
    /// Prints the vector to the console.
    /// </summary>
    public void Print(string s = "")
    {
        Console.Write(s);
        Console.WriteLine(ToString());
    }
}
