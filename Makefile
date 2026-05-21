.PHONY: main par test clean
FLAGS = -std=c++20 -pthread -Wall -Iinclude
CC = g++
MPICC = mpic++

main_seq:
	@$(CC) src/Lattice.cpp src/main_seq.cpp -o main_seq $(FLAGS) && ./main_seq

main_par: 
	@$(MPICC) -DUSE_MPI src/Lattice.cpp src/main_par.cpp -o main_par $(FLAGS)

test:
	@$(CC) src/Lattice.cpp test/test_lattice.cpp -o test_bin $(FLAGS) && ./test_bin

clean:
	@rm -f main_seq main_par test_bin
