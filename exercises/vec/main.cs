using static System.Console;
using static System.Math;

static class main
{
    // Extension method to print a double value with an optional label
    public static void print(this double x, string s = "") { Write(s); WriteLine(x); }

    // Extension method to print a Vec object with an optional label
    public static void print(this Vec v, string s = "")
    {
        Write(s);
        WriteLine($"{v.x} {v.y} {v.z}");
    }

    static int Main()
    {
        var rnd = new System.Random();

        // Generate two random vectors
        var u = new Vec(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble());
        var v = new Vec(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble());

        // Print the generated vectors
        u.print("u=");
        v.print("v=");

        WriteLine($"u={u}");
        WriteLine($"v={v}");
        WriteLine();

        Vec t;

        // Test unary negation (-u)
        t = new Vec(-u.x, -u.y, -u.z);
        (-u).print("-u =");
        t.print("t  =");
        if (Vec.Approx(t, -u)) WriteLine("test 'unary -' passed\n");

        // Test vector subtraction (u - v)
        t = new Vec(u.x - v.x, u.y - v.y, u.z - v.z);
        (u - v).print("u-v =");
        t.print("t   =");
        if (Vec.Approx(t, u - v)) WriteLine("test 'operator-' passed\n");

        // Test vector addition (u + v)
        t = new Vec(u.x + v.x, u.y + v.y, u.z + v.z);
        (u + v).print("u+v =");
        t.print("t   =");
        if (Vec.Approx(t, u + v)) WriteLine("test 'operator+' passed\n");

        // Test scalar multiplication (u * c)
        double c = rnd.NextDouble();
        t = new Vec(u.x * c, u.y * c, u.z * c);
        var tmp = u * c; // Workaround for bug in mcs
        tmp.print("u*c =");
        t.print("t   =");
        if (Vec.Approx(t, u * c)) WriteLine("test 'operator*' passed\n");

        // Test dot product (u . v)
        double d = u.x * v.x + u.y * v.y + u.z * v.z;
        d.print("u%v=");
        d.print("d  =");
        if (Vec.Approx(d, u.Dot(v))) WriteLine("test 'dot product' passed\n");

        return 0;
    }
}