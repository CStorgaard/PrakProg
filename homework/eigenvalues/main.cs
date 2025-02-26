using System;

class Program
{
    static void Main()
    {
        Console.WriteLine("Exercise A");
        Console.WriteLine("Checking Jacobi diagonalization");

        // Create a symmetric matrix (note: for eigenvalue decomposition the matrix must be symmetric)
        Mat symMatrix = new Mat(2, 2);
        symMatrix.SetElement(0, 0, 1.0);
        symMatrix.SetElement(0, 1, 2.0);
        symMatrix.SetElement(1, 0, 2.0);  // Ensure symmetry: element (1,0) equals element (0,1)
        symMatrix.SetElement(1, 1, 4.0);

        Console.WriteLine("Input Matrix:");
        symMatrix.Print();

        // Use the cyclic method to perform Jacobi eigenvalue diagonalization
        (Vec eigenvalues, Mat eigenvectors) = Mat.cyclic(symMatrix);

        Console.WriteLine("Eigenvalues:");
        eigenvalues.Print();

        Console.WriteLine("Normalized eigenvectors:");
        eigenvectors.Print();

        /// Check that V^T * A * V = D	
        Mat diagMatrix = eigenvectors.Transpose() * symMatrix * eigenvectors;
        Console.WriteLine("Diagonal matrix:");
        diagMatrix.Print();

        /// Check that V*D*V^T = A
        Mat originalMatrix = eigenvectors * diagMatrix * eigenvectors.Transpose();
        Console.WriteLine("Reconstructed matrix:");
        originalMatrix.Print();

        /// Check that V*V^T = I and V^T*V = I
        Mat identityMatrix = eigenvectors * eigenvectors.Transpose();
        Console.WriteLine("V*V^T:");
        identityMatrix.Print();

        identityMatrix = eigenvectors.Transpose() * eigenvectors;
        Console.WriteLine("V^T*V:");
        identityMatrix.Print();
    }
}
