# H-Bonding-networks
Hydrogen bonding networks from electronegativity atom (X=O, N, S) to other electronegativity atoms via water (solvent) bridge.
Use C++ program to find the HB-network between desired X atoms.
It will provide you the number of water (solvent) molecules between X1 to X2.

X1 -------- X2
X1--(H2O)n----X2



g++ -Wall -O3 -o proton_path.x proton_path_subu.cpp -std=c++11

 ./proton_path_ellipticine.x POSITION.xyz output 0 1

 The last two values are starting and ending nodes index: "0" is NH and "1" is N.
