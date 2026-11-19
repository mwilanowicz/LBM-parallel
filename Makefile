.PHONY: cavity_seq cavity_par karman_seq karman_par test clean
FLAGS = -std=c++20 -pthread -Wall -Iinclude -Wno-unused-command-line-argument -Wl,-w
CC = g++
MPICC = mpic++

# --- Cavity ---
cavity_seq:
	@$(CC) src/Lattice.cpp src/cavity/main_seq.cpp -o cavity_seq $(FLAGS) && ./cavity_seq

cavity_par: 
	@$(MPICC) -DUSE_MPI src/Lattice.cpp src/cavity/main_par.cpp -o cavity_par $(FLAGS)

# --- Karman ---
karman_seq:
	@$(CC) src/Lattice.cpp src/karman/main_seq.cpp -o karman_seq $(FLAGS) && ./karman_seq

karman_par: 
	@$(MPICC) -DUSE_MPI src/Lattice.cpp src/karman/main_par.cpp -o karman_par $(FLAGS)

# --- Test & Clean ---
test:
	@$(CC) src/Lattice.cpp test/test_lattice.cpp -o test_bin $(FLAGS) && ./test_bin

clean:
	@rm -f cavity_seq cavity_par karman_seq karman_par test_bin
