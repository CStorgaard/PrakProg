using System;
using System.Text;

/// <summary>
/// Represents a matrix and provides various matrix operations.
/// </summary>
public class Mat
{
    private double[,] data;

    /// <summary>
    /// Number of rows in the matrix.
    /// </summary>
    public int Rows { get; }

    /// <summary>
    /// Number of columns in the matrix.
    /// </summary>
    public int Cols { get; }

    /// <summary>
    /// Constructs a matrix with the specified number of rows and columns.
    /// All elements are initialized to zero.
    /// </summary>
    public Mat(int rows, int cols)
    {
        if (rows <= 0 || cols <= 0)
            throw new ArgumentException("Matrix dimensions must be positive integers.");
        Rows = rows;
        Cols = cols;
        data = new double[rows, cols];
    }

    /// <summary>
    /// Constructs a matrix using a 2D array.
    /// </summary>
    public Mat(double[,] array)
    {
        if (array == null)
            throw new ArgumentNullException(nameof(array));
        Rows = array.GetLength(0);
        Cols = array.GetLength(1);
        data = (double[,])array.Clone();
    }

    /// <summary>
    /// Gets or sets the matrix element at the specified row and column.
    /// </summary>
    public double this[int row, int col]
    {
        get { return data[row, col]; }
        set { data[row, col] = value; }
    }

    /// <summary>
    /// Sets the element at the specified row and column.
    /// </summary>
    public void SetElement(int row, int col, double value)
    {
        this[row, col] = value;
    }

    /// <summary>
    /// Adds two matrices of the same dimensions.
    /// </summary>
    public static Mat operator +(Mat a, Mat b)
    {
        if (a.Rows != b.Rows || a.Cols != b.Cols)
            throw new ArgumentException("Matrix dimensions must match for addition.");
        Mat result = new Mat(a.Rows, a.Cols);
        for (int i = 0; i < a.Rows; i++)
            for (int j = 0; j < a.Cols; j++)
                result[i, j] = a[i, j] + b[i, j];
        return result;
    }

    /// <summary>
    /// Subtracts matrix b from matrix a.
    /// </summary>
    public static Mat operator -(Mat a, Mat b)
    {
        if (a.Rows != b.Rows || a.Cols != b.Cols)
            throw new ArgumentException("Matrix dimensions must match for subtraction.");
        Mat result = new Mat(a.Rows, a.Cols);
        for (int i = 0; i < a.Rows; i++)
            for (int j = 0; j < a.Cols; j++)
                result[i, j] = a[i, j] - b[i, j];
        return result;
    }

    /// <summary>
    /// Multiplies a matrix by a scalar.
    /// </summary>
    public static Mat operator *(Mat a, double scalar)
    {
        Mat result = new Mat(a.Rows, a.Cols);
        for (int i = 0; i < a.Rows; i++)
            for (int j = 0; j < a.Cols; j++)
                result[i, j] = a[i, j] * scalar;
        return result;
    }

    /// <summary>
    /// Multiplies a scalar by a matrix.
    /// </summary>
    public static Mat operator *(double scalar, Mat a) => a * scalar;

    /// <summary>
    /// Multiplies two matrices (matrix multiplication).
    /// </summary>
    public static Mat operator *(Mat a, Mat b)
    {
        if (a.Cols != b.Rows)
            throw new ArgumentException("Invalid matrix dimensions for multiplication.");
        Mat result = new Mat(a.Rows, b.Cols);
        for (int i = 0; i < a.Rows; i++)
        {
            for (int j = 0; j < b.Cols; j++)
            {
                double sum = 0;
                for (int k = 0; k < a.Cols; k++)
                {
                    sum += a[i, k] * b[k, j];
                }
                result[i, j] = sum;
            }
        }
        return result;
    }

    /// <summary>
    /// Multiplies a matrix with a vector.
    /// </summary>
    public static Vec operator *(Mat m, Vec v)
    {
        if (m.Cols != v.Length)
            throw new ArgumentException("Matrix columns must match vector length.");
        Vec result = new Vec(m.Rows);
        for (int i = 0; i < m.Rows; i++)
        {
            double sum = 0;
            for (int j = 0; j < m.Cols; j++)
            {
                sum += m[i, j] * v[j];
            }
            result[i] = sum;
        }
        return result;
    }

    /// <summary>
    /// Returns the negation of the matrix.
    /// </summary>
    public static Mat operator -(Mat a) => a * (-1);

    /// <summary>
    /// Returns the transpose of the matrix.
    /// </summary>
    public Mat Transpose()
    {
        Mat result = new Mat(Cols, Rows);
        for (int i = 0; i < Rows; i++)
            for (int j = 0; j < Cols; j++)
                result[j, i] = this[i, j];
        return result;
    }

    /// <summary>
    /// Get a column of the matrix as a vector.
    /// </summary>
    public Vec GetCol(int col)
    {
        if (col < 0 || col >= Cols)
            throw new ArgumentOutOfRangeException(nameof(col));
        Vec v = new Vec(Rows);
        for (int i = 0; i < Rows; i++)
            v[i] = this[i, col];
        return v;
    }

    /// <summary>
    /// Get a row of the matrix as a vector.
    /// </summary>
    public Vec GetRow(int row)
    {
        if (row < 0 || row >= Rows)
            throw new ArgumentOutOfRangeException(nameof(row));
        Vec v = new Vec(Cols);
        for (int j = 0; j < Cols; j++)
            v[j] = this[row, j];
        return v;
    }

    /// <summary>
    /// Set the column to the values of the vector.
    /// </summary>
    public void SetCol(int col, Vec v)
    {
        if (col < 0 || col >= Cols)
            throw new ArgumentOutOfRangeException(nameof(col));
        if (v.Length != Rows)
            throw new ArgumentException("Vector length must match matrix row count.");
        for (int i = 0; i < Rows; i++)
            this[i, col] = v[i];
    }

    /// <summary>
    /// Perform QR decomposition of the matrix using the classical Gram–Schmidt process.
    /// This method uses out parameters.
    /// </summary>
    public void QRDecomposition(out Mat Q, out Mat R)
    {
        Q = new Mat(Rows, Cols);
        R = new Mat(Cols, Cols);
        // Create a copy of the current matrix.
        Mat A = new Mat(data);
        for (int j = 0; j < Cols; j++)
        {
            Vec v = A.GetCol(j);
            for (int i = 0; i < j; i++)
            {
                Vec q = Q.GetCol(i);
                R[i, j] = q.Dot(v);
                v = v - q * R[i, j];
            }
            R[j, j] = v.Norm();
            Q.SetCol(j, v / R[j, j]);
        }
    }

    /// <summary>
    /// Overload of QRDecomposition that returns a tuple containing Q and R.
    /// </summary>
    public (Mat Q, Mat R) QRDecomposition()
    {
        Mat Q, R;
        QRDecomposition(out Q, out R);
        return (Q, R);
    }

    /// <summary>
    /// Solves the system of linear equations Ax = b using QR decomposition.
    /// </summary>
    public Vec SolveQR(Vec b)
    {
        if (Rows != Cols || Rows != b.Length)
            throw new ArgumentException("Matrix dimensions must match vector length.");
        Mat Q, R;
        QRDecomposition(out Q, out R);
        Vec y = Q.Transpose() * b;
        Vec x = new Vec(Cols);
        for (int i = Cols - 1; i >= 0; i--)
        {
            double sum = 0;
            for (int j = i + 1; j < Cols; j++)
                sum += R[i, j] * x[j];
            x[i] = (y[i] - sum) / R[i, i];
        }
        return x;
    }

    /// <summary>
    /// Computes the determinant of the matrix using QR decomposition.
    /// </summary>
    public double Determinant()
    {
        if (Rows != Cols)
            throw new InvalidOperationException("Matrix must be square.");
        Mat Q, R;
        QRDecomposition(out Q, out R);
        double det = 1;
        for (int i = 0; i < Cols; i++)
            det *= R[i, i];
        return det;
    }

    /// <summary>
    /// Computes the inverse of the matrix using QR decomposition.
    /// </summary>
    public Mat Inverse()
    {
        if (Rows != Cols)
            throw new InvalidOperationException("Matrix must be square.");
        Mat Q, R;
        QRDecomposition(out Q, out R);
        Mat inv = new Mat(Cols, Cols);
        // Solve Ax = e_i for each column of the identity matrix.
        for (int i = 0; i < Cols; i++)
        {
            Vec b = new Vec(Cols);
            b[i] = 1;
            Vec x = SolveQR(b);
            for (int j = 0; j < Cols; j++)
                inv[j, i] = x[j];
        }
        return inv;
    }


    /// <summary>
    /// Multiplies matrix A on the right by the Jacobi rotation matrix J.
    /// That is, it updates columns p and q of A.
    /// </summary>
    public static void timesJ(Mat A, int p, int q, double theta)
    {
        double c = Math.Cos(theta);
        double s = Math.Sin(theta);
        for (int i = 0; i < A.Rows; i++)
        {
            double temp = A[i, p];
            double temp2 = A[i, q];
            A[i, p] = c * temp - s * temp2;
            A[i, q] = s * temp + c * temp2;
        }
    }
    
    /// <summary>
    /// Multiplies matrix A on the left by the transpose of the Jacobi rotation matrix, J^T.
    /// That is, it updates rows p and q of A.
    /// </summary>
    public static void Jtimes(Mat A, int p, int q, double theta)
    {
        double c = Math.Cos(theta);
        double s = Math.Sin(theta);
        for (int j = 0; j < A.Cols; j++)
        {
            double temp = A[p, j];
            double temp2 = A[q, j];
            A[p, j] = c * temp + s * temp2;
            A[q, j] = -s * temp + c * temp2;
        }
    }
    
    /// <summary>
    /// Returns a copy of the matrix A.
    /// </summary>
    public static Mat Copy(Mat A)
    {
        Mat copy = new Mat(A.Rows, A.Cols);
        for (int i = 0; i < A.Rows; i++)
        {
            for (int j = 0; j < A.Cols; j++)
            {
                copy[i, j] = A[i, j];
            }
        }
        return copy;
    }

    /// <summary>
    /// Returns the identity matrix of size n.
    /// </summary>
    public static Mat Identity(int n)
    {
        Mat I = new Mat(n, n);
        for (int i = 0; i < n; i++)
        {
            I[i, i] = 1;
        }
        return I;
    }

    

    /// <summary>
    /// Performs cyclic Jacobi diagonalization of a symmetric matrix.
    /// </summary>
    public static (Vec, Mat) cyclic(Mat M)
    {
        int n = M.Rows;
        // Work on a copy so as not to modify the original matrix.
        Mat A = Copy(M);
        // Initialize V to the identity matrix.
        Mat V = Identity(n);
        
        double eps = 1e-9;
        int maxIterations = 100 * n * n;
        int iter = 0;
        bool changed = true;
        
        while (changed && iter < maxIterations)
        {
            changed = false;
            for (int p = 0; p < n - 1; p++)
            {
                for (int q = p + 1; q < n; q++)
                {
                    double Apq = A[p, q];
                    if (Math.Abs(Apq) > eps)
                    {
                        double App = A[p, p];
                        double Aqq = A[q, q];
                        double tau = (Aqq - App) / (2 * Apq);
                        double t = Math.Sign(tau) / (Math.Abs(tau) + Math.Sqrt(1 + tau * tau));
                        double c = 1.0 / Math.Sqrt(1 + t * t);
                        double s = t * c;
                        
                        // Update diagonal elements.
                        double App_new = c * c * App - 2 * s * c * Apq + s * s * Aqq;
                        double Aqq_new = s * s * App + 2 * s * c * Apq + c * c * Aqq;
                        A[p, p] = App_new;
                        A[q, q] = Aqq_new;
                        A[p, q] = 0;
                        A[q, p] = 0;
                        
                        // Update remaining elements.
                        for (int r = 0; r < n; r++)
                        {
                            if (r != p && r != q)
                            {
                                double Arp = A[r, p];
                                double Arq = A[r, q];
                                double newArp = c * Arp - s * Arq;
                                double newArq = s * Arp + c * Arq;
                                A[r, p] = newArp;
                                A[p, r] = newArp; // preserve symmetry
                                A[r, q] = newArq;
                                A[q, r] = newArq; // preserve symmetry
                            }
                        }
                        
                        // Update the eigenvector matrix V.
                        for (int r = 0; r < n; r++)
                        {
                            double Vrp = V[r, p];
                            double Vrq = V[r, q];
                            double newVrp = c * Vrp - s * Vrq;
                            double newVrq = s * Vrp + c * Vrq;
                            V[r, p] = newVrp;
                            V[r, q] = newVrq;
                        }
                        
                        changed = true;
                    }
                }
            }
            iter++;
        }
        
        // Extract the eigenvalues from the diagonal of A.
        Vec w = new Vec(n);
        for (int i = 0; i < n; i++)
        {
            w[i] = A[i, i];
        }
        
        return (w, V);
    }


    /// <summary>
    /// Checks if two matrices are approximately equal element-wise.
    /// </summary>
    public bool Approx(Mat other, double acc = 1e-9, double eps = 1e-9)
    {
        if (other == null || this.Rows != other.Rows || this.Cols != other.Cols)
            return false;
        for (int i = 0; i < Rows; i++)
        {
            for (int j = 0; j < Cols; j++)
            {
                if (!ApproxEqual(this[i, j], other[i, j], acc, eps))
                    return false;
            }
        }
        return true;
    }

    /// <summary>
    /// Helper method to compare two doubles for approximate equality.
    /// </summary>
    public static bool ApproxEqual(double a, double b, double acc = 1e-9, double eps = 1e-9)
    {
        if (Math.Abs(a - b) < acc)
            return true;
        if (Math.Abs(a - b) < (Math.Abs(a) + Math.Abs(b)) * eps)
            return true;
        return false;
    }

    /// <summary>
    /// Returns a string representation of the matrix.
    /// </summary>
    public override string ToString()
    {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < Rows; i++)
        {
            for (int j = 0; j < Cols; j++)
            {
                sb.Append(this[i, j].ToString("G"));
                if (j < Cols - 1)
                    sb.Append(" ");
            }
            if (i < Rows - 1)
                sb.AppendLine();
        }
        return sb.ToString();
    }

    /// <summary>
    /// Prints the matrix to the console.
    /// </summary>
    public void Print(string s = "")
    {
        Console.Write(s);
        Console.WriteLine(ToString());
    }
}
