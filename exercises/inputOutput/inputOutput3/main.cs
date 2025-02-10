class main{
    static int Main(string[] args){
        foreach(string arg in args) System.Console.Out.WriteLine($"sdtout: {arg}");
        foreach(string arg in args) System.Console.Error.WriteLine($"stderr: {arg}");
        string inputfile = " ", outputfile = " ";
        foreach(string arg in args){
            string[] words = arg.Split(":");
            if(words[0] == "-input")inputfile = words[1];
            if(words[0] == "-output")outputfile = words[1];
        }

        System.Console.Error.WriteLine($"inputfile: {inputfile} outputfile: {outputfile}");

    if(inputfile == " " || outputfile == " ") return 0;

        var instream = new System.IO.StreamReader(inputfile);
        var outstream = new System.IO.StreamWriter(outputfile, append: true);
        
        for(
            string line = instream.ReadLine();
            line != null;
            line = instream.ReadLine()
        ){
            outstream.WriteLine($"line for instream: {line}");
        }

        instream.Close();
        outstream.Close();
        return 0;
    }
}