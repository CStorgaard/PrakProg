.PHONY: all run-tests test clean

# Default target builds main.exe and runs the tests.
all: run-tests

# run-tests builds Out.txt by running several tests.
Out.txt: Makefile main.exe input.txt
	@echo "Truncating log file..."
	@> log
	@echo "hello" > $@
	@echo "1,2 3,4 5,6" | sed --expression="s/,/\./g" >> $@
	@echo "hello world from stdin" | mono main.exe 23,23, 44,44 55,55 1>> $@ 2> log
	@cat input.txt | mono main.exe 1>> $@ 2>> log
	@mono main.exe < input.txt 1>> $@ 2>> log
	@echo "" | mono main.exe -input:input.txt -output:$@ 2>> log

# Generate input.txt from the Makefile
input.txt: Makefile
	@echo "these go $@" > $@

# Compile the main program
main.exe: main.cs
	mcs -target:exe -out:$@ $<

# An example test target (you can add more commands as needed)
test:
	@echo "Test target: \$@"
	@echo '$@'
	@echo "$@"

# Clean up generated files
clean:
	rm -f *.exe *.txt log

# Convenience target for running tests
run-tests: Out.txt
