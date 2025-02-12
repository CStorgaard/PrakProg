using System;
using static System.Console;

public class genlist<T>
{
    public T[] data; // Field to store the items

    // Property to return the current number of items.
    public int size => data.Length;

    // Indexer to allow access via [index]
    public T this[int i] => data[i];

    // Constructor initializes an empty array.
    public genlist() {
        data = new T[0];
    }

    // Method to add an item by allocating a new array, copying the old data,
    // and then adding the new item.
    public void add(T item) {
        T[] newdata = new T[size + 1];
        Array.Copy(data, newdata, size);
        newdata[size] = item;
        data = newdata;
    }

    // Method to remove the element at index i.
    public void remove(int i) {
        if (i < 0 || i >= size) {
            throw new ArgumentOutOfRangeException(nameof(i), "Index out of range.");
        }
        T[] newdata = new T[size - 1];
        // Copy elements before index i.
        if (i > 0)
            Array.Copy(data, 0, newdata, 0, i);
        // Copy elements after index i.
        if (i < size - 1)
            Array.Copy(data, i + 1, newdata, i, size - i - 1);
        data = newdata;
    }
}

public class Program {
    public static void Main(string[] args) {
        // Create an instance of the generic list to hold arrays of doubles.
        var list = new genlist<double[]>();

        // Delimiters and options for splitting input lines.
        char[] delimiters = { ' ', '\t' };
        var options = StringSplitOptions.RemoveEmptyEntries;

        // Read input lines until null is encountered.
        for (string line = ReadLine(); line != null; line = ReadLine()) {
            var words = line.Split(delimiters, options);
            int n = words.Length;
            var numbers = new double[n];
            for (int i = 0; i < n; i++) {
                numbers[i] = double.Parse(words[i]);
            }
            list.add(numbers);
        }

        // For example, remove the second element (index 1) if it exists.
        if (list.size > 1) {
            list.remove(1);
        }

        // Output each array of numbers.
        for (int i = 0; i < list.size; i++) {
            var numbers = list[i];
            foreach (var number in numbers) {
                Write($"{number:0.00e+00;-0.00e+00} ");
            }
            WriteLine();
        }
    }
}
