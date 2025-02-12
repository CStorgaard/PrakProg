.PHONY: all clean

# The default target builds main.exe and then runs the program to create Out.txt
all: Out.txt

# Out.txt depends on main.exe and input.txt
Out.txt: main.exe input.txt
	@echo "Running main.exe with input from input.txt..."
	@echo "Removing third column from input and writing to Out.txt..."
	@mono main.exe < input.txt > Out.txt

# Compile main.exe from main.cs
main.exe: main.cs
	mcs -target:exe -out:main.exe main.cs

# Create an input file containing the table with tabs and spaces
input.txt:
	@printf "129.24 24.8             4.847\n\t\t88.6   33.745 7.78\n30.39  99.811              6.723\n\t\t-1.33   96.3   2.23\n" > input.txt

# Clean up generated files
clean:
	rm -f main.exe Out.txt input.txt
