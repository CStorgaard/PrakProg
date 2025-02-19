TARGET = main.exe
SOURCE = main.cs
N = 1e8
TIME = time --portability --append --output $@

# Default target builds main.exe and runs all tests.
all: $(TARGET) Out.txt Out_parallel.txt Out_threadlocal.txt

# Compile main.cs into main.exe using mcs.
$(TARGET): $(SOURCE)
	mcs -out:$(TARGET) $(SOURCE)

# Manual threading version output with multiple thread counts.
Out.txt: $(TARGET) Makefile
	@echo "Manual threading:" > $@
	@echo "Using 1 thread:" >> $@
	$(TIME) mono $< -terms:$(N) -threads:1 >> $@
	@echo "Using 2 threads:" >> $@
	$(TIME) mono $< -terms:$(N) -threads:2 >> $@
	@echo "Using 3 threads:" >> $@
	$(TIME) mono $< -terms:$(N) -threads:3 >> $@
	@echo "Using 4 threads:" >> $@
	$(TIME) mono $< -terms:$(N) -threads:4 >> $@

# Flawed Parallel.For version (shared variable race condition).
Out_parallel.txt: $(TARGET) Makefile
	@echo "Flawed Parallel.For version (shared variable):" > $@
	$(TIME) mono $< -terms:$(N) -parallel >> $@

# Parallel.For version using ThreadLocal.
Out_threadlocal.txt: $(TARGET) Makefile
	@echo "Parallel.For with ThreadLocal:" > $@
	$(TIME) mono $< -terms:$(N) -threadlocal >> $@

# Clean up generated files.
clean:
	rm -f $(TARGET) Out.txt Out_parallel.txt Out_threadlocal.txt
