using static System.Math;
using System;

class Program
{
    static void Main()
    {
        Console.WriteLine("Exercise A");
        Console.WriteLine("Checking if decomp works for tall matrix");

        // Use "Mat" (the library class) instead of "Matrix"
        Mat tallMatrix = new Mat(4, 2); // Creating a tall matrix with 4 rows and 2 columns

        // Setting elements of the matrix using SetElement
        tallMatrix.SetElement(0, 0, 1.0);
        tallMatrix.SetElement(0, 1, 2.0);
        tallMatrix.SetElement(1, 0, 3.0);
        tallMatrix.SetElement(1, 1, 4.0);
        tallMatrix.SetElement(2, 0, 5.0);
        tallMatrix.SetElement(2, 1, 6.0);
        tallMatrix.SetElement(3, 0, 7.0);
        tallMatrix.SetElement(3, 1, 8.0);

        // Printing the matrix
        tallMatrix.Print();

        // Factorize it using QR decomposition (using the new overload)
        var qr = tallMatrix.QRDecomposition();

        // Printing the Q matrix and the R matrix
        qr.Q.Print("Q = ");
        qr.R.Print("R = ");

        // Check that Q^T * Q = I
        var qtq = qr.Q.Transpose() * qr.Q;

        // Check that Q*R = A
        var qrCheck = qr.Q * qr.R;

        // Print the results
        qtq.Print("Q^T * Q = ");
        qrCheck.Print("Q * R = ");


        Console.WriteLine("Checking if solve works on square matrix");

        // Create a square matrix
        Mat squareMatrix = new Mat(2, 2);

        // Set elements of the matrix
        squareMatrix.SetElement(0, 0, 1.0);
        squareMatrix.SetElement(0, 1, 2.0);
        squareMatrix.SetElement(1, 0, 3.0);
        squareMatrix.SetElement(1, 1, 4.0);

        // Print the matrix
        squareMatrix.Print("A = ");

        // Create a vector
        Vec b = new Vec(2);

        // Set elements of the vector
        b.SetElement(0, 5.0);
        b.SetElement(1, 6.0);

        // Print the vector
        b.Print("b = ");

        // Solve the system
        Vec x = squareMatrix.SolveQR(b);

        // Print the solution
        x.Print("x = ");

        // Check that the solution is correct
        Vec bCheck = squareMatrix * x;
        bCheck.Print("A * x = "); // Should be equal to b

        Console.WriteLine("Exercise B");

        // Create a square matrix
        Mat squareMatrix2 = new Mat(2, 2);

        // Set elements of the matrix
        squareMatrix2.SetElement(0, 0, 1.0);
        squareMatrix2.SetElement(0, 1, 2.0);
        squareMatrix2.SetElement(1, 0, 3.0);
        squareMatrix2.SetElement(1, 1, 4.0);

        // Print the matrix
        squareMatrix2.Print("A = ");

        // Factorize it using QR decomposition
        var qr2 = squareMatrix2.QRDecomposition();

        // Print the Q matrix and the R matrix
        qr2.Q.Print("Q = ");
        qr2.R.Print("R = ");

        // Calculate the inverse of the matrix A
        Mat invA = squareMatrix2.Inverse();

        // Print the inverse
        invA.Print("A^-1 = ");

        // Check that A * A^-1 = I

        // Calculate the product
        Mat check = squareMatrix2 * invA;

        // Print the result
        check.Print("A * A^-1 = ");
    }
}
