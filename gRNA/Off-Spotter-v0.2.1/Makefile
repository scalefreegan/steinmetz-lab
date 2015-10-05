#Off-Spotter make file
all: Results Load_Memory Table_Creation Detach_Memory

Results: Off-Spotter_Results.cpp
	g++ Off-Spotter_Results.cpp -m64 -O3 -o Results

Load_Memory: Off-Spotter_Load_Memory.cpp
	g++ Off-Spotter_Load_Memory.cpp -m64 -O3 -o Load_Memory

Table_Creation: Off-Spotter_Table_Creation.cpp
	g++ Off-Spotter_Table_Creation.cpp -m64 -O3 -o Table_Creation

Detach_Memory: Off-Spotter_Detach_Memory.cpp
	g++ Off-Spotter_Detach_Memory.cpp -m64 -O3 -o Detach_Memory

clean:
	rm Detach_Memory Table_Creation Load_Memory Results