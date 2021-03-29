// basic file operations
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <time.h>
#include <float.h>
#include <set>
#include <sstream>

using namespace std;
const int DIM = 3;
const int INF = 10000;
// int *through;
// int *dist;
const double pi=atan(1.0)*4.0;
double dist_calc(double L[][DIM], double *x, double *y);
double angle_calc(double L[][DIM], double *x, double *y, double *z);
double scalar_product(double *x, double *y);
double norm(double *x);
void printUsage(string program);

int main (int argc, char * argv[]) {
    //Defining time variables to examine the total computing time
    double t_start, t_stop;
    t_start=clock();

    //Reading the command lines
    
	if (argc != 5) {printUsage(argv[0]); return -1;};

    string fn_xyz=argv[1];
    string fn_out=argv[2];
    string s_str=argv[3];
    string g_str=argv[4];
    int s = stoi(s_str);
    int g = stoi(g_str);
    cout << " >>>>>>>>>>> input parameters " << endl;
    cout << "file name to be read                 : " << fn_xyz << endl;
    cout << "output file for ...                  : " << fn_out << endl;
    cout << "starting node of proton wire: only for N: " << s << endl;
    cout << "ending node of proton wire: only for N: " << g << endl;
    cout << endl;
    
    //Total time fames to be used for the calcuation
    // int tot_frames = (t_f-t_i)/dt+1;
    cout << " >>>>>>>>>>> xyz file format!" << endl;
    static const int n_head = 9;     //# of header of file to be read before the coordinates
    cout << " Number of words in second line : " << n_head << endl;
    //Number of atoms of protein and water at each frame
    static const int n_atoms = 294;      //number of the total atoms in each frame for neutral system
    static const int n_pro = 33;		    //# of protein Urea atoms
    static const int n_pro_h = 1; 		//# of H of NH of Urea
    static const int n_pro_n = 2;		//# of N of NH of Urea
    static const int n_pro_o = 0; 		//# of O of CO of Urea
    static const int n_ow = 87;     		//# of OW atoms
    static const int n_hw = 174;   		//# of HW atoms for the system
    int pro_n_ndx[n_pro_n] = {6, 6}; //index of N of Urea starting from 0
    int pro_h_ndx[n_pro_h] = {7};  //index of O of Urea starting from 0
    int pro_o_ndx[n_pro_o] = {};  //index of O of Urea starting from 0
    cout << " Total number of atoms in the system : " << n_atoms << endl;
    cout << " Total number of atoms of protein    : " << n_pro << endl;
    cout << " Total number of nitrogen of NH group: " << n_pro_n << endl;
    for (int i=0; i<n_pro_n; ++i){ cout << i << " th index of N of NH >> " << pro_n_ndx[i] << endl;}
    cout << " Total number of hydrogen of NH group: " << n_pro_h << endl;
    for (int i=0; i<n_pro_h; ++i){ cout << i << " th index of H of NH >> " << pro_h_ndx[i] << endl;}
    cout << " Total number of oxygen of protein : " << n_pro_o << endl;
    for (int i=0; i<n_pro_o; ++i){ cout << i << " th index of O of CO >> " << pro_o_ndx[i] << endl;}
    cout << " Total number of water oxygens      : " << n_ow << endl;
    cout << " Total number of water hydroggens   : " << n_hw << endl;
    // Dynamically allocate coordinates arrays
    double **pro_pos, **ow_pos, **h_pos; // array of the coordinates at each frame
    pro_pos = new double *[n_pro]();  //protein coordinates
    ow_pos = new double *[n_ow]();    //OW coordinates
    h_pos = new double *[n_hw+n_pro_h]();  //HW coordinates and H of NH (so we plus 4)
    // Dynamically allocate index of OW and HW at each frame
    int *ow_ndx, *h_ndx;
    ow_ndx = new int [n_ow]();    //OW coordinates
    h_ndx = new int [n_hw+n_pro_h]();  //HW coordinates and H of NH (so we plus 4)
    
    for (int i=0; i<n_pro; i++){
        pro_pos[i]= new double[DIM]();
    }
    for (int i=0; i<n_ow; i++){
        ow_pos[i]= new double[DIM]();
    }
    for (int i=0; i<n_hw+n_pro_h; i++){
        h_pos[i]= new double[DIM]();
    }
	
    string fn_temp;
    //reading input files "---.xyz"
    ifstream traj_file(fn_xyz.c_str());
    //writing output file "---.out" to print sum data
    fn_temp = fn_out + ".sum";
    ofstream sum_file(fn_temp.c_str());
    fn_temp = fn_out + ".path";
    ofstream path_file(fn_temp.c_str());
    fn_temp = fn_out + ".length";
    ofstream length_file(fn_temp.c_str());
    fn_temp = fn_out + ".in";
    ofstream in_file(fn_temp.c_str());
    fn_temp = fn_out + ".out";
    ofstream out_file(fn_temp.c_str());
    // fn_temp = fn_out + ".h_chk";
    // ofstream proton_chk_file(fn_temp.c_str());
    
    if(!traj_file) {
        cout << "Cannot open input file." << endl;
        return 1;
    }
	
    int n_atoms1, n_count_ow, n_count_hw, n_count_atom;
    string atom_name;
    string str;
    int n_frame; //frame index

    n_frame = 1;
    while (traj_file.is_open()){
        // for (int i=0; i<tot_frames; ++i){
        /*reading the header files of each frame*/
        n_atoms1 = 0;
        traj_file >> n_atoms1;
        if((n_atoms1 != n_atoms) && (n_frame == 0)) {
            cout << "total number of atoms is different from the given data! pls check!" << endl;
            exit (EXIT_FAILURE);
        }else if((n_atoms1 != n_atoms) && (n_frame > 0)) {
            cout << "End of Calculation! Good Job!" << endl;
            break;
        }else{
            // cout << "# of frame    =  "<< n_frame << endl;
            // cout << "number of atoms:" << n_atoms1 << endl;
        }
        for (int i=0; i<n_head; ++i){
            traj_file >> str;
             // cout << "second line " << i << ": " << str << endl;
        }
        /*reading the protein coordinates of each frame*/
        for(int j=0; j<n_pro; ++j){
            traj_file >> atom_name;
            for(int k=0; k<DIM; ++k){
                traj_file >> pro_pos[j][k];
            }
            // cout << atom_name << " " << j << " " << pro_pos[j][0] << " " << pro_pos[j][1] << " " << pro_pos[j][2] << endl;
        }

        n_count_ow = 0; //Set starting index of O atom as 0
        n_count_hw = 0; //Set starting index of H atom as 0
        n_count_atom = n_pro;
        for(int j=0; j<(n_ow+n_hw); ++j){
            traj_file >> atom_name;
            if (atom_name == "O"){
                for(int k=0; k<DIM; ++k){
                    traj_file >> ow_pos[n_count_ow][k];
                }
                ow_ndx[n_count_ow] = n_count_atom;
                // cout << atom_name <<  " # " << n_count_atom << " " << n_count_ow << ":" << ow_pos[n_count_ow][0] << " " << ow_pos[n_count_ow][1] << " " << ow_pos[n_count_ow][2] << endl;
                n_count_ow += 1;
				
            }else if(atom_name == "H"){
                for(int k=0; k<DIM; ++k){
                    traj_file >> h_pos[n_count_hw][k];
                }
                h_ndx[n_count_hw] = n_count_atom;
                // cout << atom_name << " # " << n_count_atom << " " << n_count_hw << ":" << h_pos[n_count_hw][0] << " " << h_pos[n_count_hw][1] << " " << h_pos[n_count_hw][2] << endl;
                n_count_hw += 1;
            }else{
                cout << "ERROR : Please check of atom name carefully 1" << endl;
                return -1;
            }
            n_count_atom += 1;
        }
		for(int j=0; j<n_pro_h; ++j){
			for(int k=0;k<DIM;++k){
				h_pos[n_hw+j][k] = pro_pos[pro_h_ndx[j]][k];  // H atom of NH group is now involved in the hydrogen atoms of water
			}
			h_ndx[n_hw+j] = pro_h_ndx[j];
			// cout << " H of NH" << pro_h_ndx[j] << " " << n_hw+j << ":" << h_pos[n_hw+j][0] << " " << h_pos[n_hw+j][1] << " " << h_pos[n_hw+j][2] << endl;
        }

        double L[DIM][DIM]={{14.861, 0, 0},
                            {0, 14.861, 0},
                            {0, 0, 14.861}};

        /* adjacent matrix for water and Urea */
        int **adjmatrix, **hb_h_ndx;            	 //  array of the adjacency matrix at each frame
        bool **conn_matrix;      	 	//  array of the connection matrix at each frame
        int V = n_ow+n_pro_n+n_pro_o;   //  number of all donors and acceptors
        adjmatrix = new int *[V]();     //  construction of adjcent matrix
        hb_h_ndx = new int *[V]();     //  construction of adjcent matrix
        conn_matrix = new bool *[V]();     //  construction of adjcent matrix
        for (int i=0; i<V; i++){
            adjmatrix[i]= new int[V]();
            hb_h_ndx[i]= new int[V]();
            conn_matrix[i]= new bool[V]();
        }

		int network_ndx[V];

        // ====> HB information <============

        double rcut_oo=3.5, rcut_oh=2.5, angcut_oho=130; //, r_bond=1.2;
        //cout << " rcut_oo: " << rcut_oo << endl;
        //cout << " rcut_oh: " << rcut_oh << endl;
        //cout << " angcut_oho: " << angcut_oho << endl;

        //====> Construction of Adjacent Matrix <=====
		
        for(int i=0; i<V; ++i){
            for(int j=0; j<V; ++j){
                adjmatrix[i][j] = INF;
                hb_h_ndx[i][j] = -1;
                conn_matrix[i][j] = 0;
            }
        }
        // ====> index to serial number <============
		
        for(int i = 0; i < n_pro_n; ++i){
            network_ndx[i] = pro_n_ndx[i]+1;
        }
	for(int i = 0; i < n_ow; ++i){
            network_ndx[n_pro_n+i] = ow_ndx[i]+1;
        }
	for(int i = 0; i < n_pro_o; ++i){
            network_ndx[n_pro_n+n_ow+i] = pro_o_ndx[i]+1;
	}

        double d1, d2, hb_ang;
		for(int i = 0; i < n_pro_n; ++i){
			for(int j = 0; j < n_ow; ++j){
          // Length of the box L //
                if(dist_calc(L, pro_pos[pro_n_ndx[i]], ow_pos[j]) < rcut_oo) { //distance between two oxygen atoms
                    for(int k = 0; k < n_hw+n_pro_h; ++k){
                        hb_ang = angle_calc(L, pro_pos[pro_n_ndx[i]], h_pos[k], ow_pos[j]);
                        if (hb_ang > angcut_oho){
                            d1 = dist_calc(L, pro_pos[pro_n_ndx[i]], h_pos[k]);
                            d2 = dist_calc(L, ow_pos[j], h_pos[k]);
							//if ((d2 > d1) && (d2 < rcut_oh) && (d1 < r_bond)) {
                            if ((d2 > d1) && (d2 < rcut_oh)) {
                                adjmatrix[i][n_pro_n+j] = 1;
                                conn_matrix[i][n_pro_n+j] = 1;            //  i -> j is connected
								hb_h_ndx[i][n_pro_n+j] = h_ndx[k]+1;
								// cout << network_ndx[i] << " " << network_ndx[n_pro_n+j] <<  " " << hb_h_ndx[i][n_pro_n+j] << " " << adjmatrix[i][n_pro_n+j] << " " << adjmatrix[n_pro_n+j][i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            //}else if((d2 < d1) && (d1 < rcut_oh) && (d2 < r_bond)) {
							}else if((d2 < d1) && (d1 < rcut_oh)){
                                adjmatrix[n_pro_n+j][i] = 1;
                                conn_matrix[n_pro_n+j][i] = 1;            //  j -> i is connected
								hb_h_ndx[n_pro_n+j][i] = h_ndx[k]+1;
								// cout << network_ndx[i] << " " << network_ndx[n_pro_n+j] <<  " " << hb_h_ndx[n_pro_n+j][i] << " " << adjmatrix[i][n_pro_n+j] << " " << adjmatrix[n_pro_n+j][i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            //}else if((d2 == d1) && (d1 < r_bond)){
							}else if((d2 == d1) && (d1 < rcut_oh)){
                                adjmatrix[i][n_pro_n+j] = 1;
                                adjmatrix[n_pro_n+j][i] = 1;
                                conn_matrix[i][n_pro_n+j] = 1;            //  i -> j is connected
                                conn_matrix[n_pro_n+j][i] = 1;            //  j -> i is connected
								hb_h_ndx[i][n_pro_n+j] = h_ndx[k]+1;
								hb_h_ndx[n_pro_n+j][i] = h_ndx[k]+1;
								// cout << network_ndx[i] << " " << network_ndx[n_pro_n+j] <<  " " << hb_h_ndx[i][n_pro_n+j] << " " << adjmatrix[i][n_pro_n+j] << " " << adjmatrix[n_pro_n+j][i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            }
							
                        }
                        
                    }
                }
            }
			for(int j = 0; j < n_pro_o; ++j){
				if(dist_calc(L, pro_pos[pro_n_ndx[i]], pro_pos[pro_o_ndx[j]]) < rcut_oo) { //distance between two oxygen atoms
                    for(int k = 0; k < n_hw+n_pro_h; ++k){
                        hb_ang = angle_calc(L, pro_pos[pro_n_ndx[i]], h_pos[k], pro_pos[pro_o_ndx[j]]);
                        if (hb_ang > angcut_oho){
                            d1 = dist_calc(L, pro_pos[pro_n_ndx[i]], h_pos[k]);
                            d2 = dist_calc(L, pro_pos[pro_o_ndx[j]], h_pos[k]);
							//if ((d2 > d1) && (d2 < rcut_oh) && (d1 < r_bond)) {
                            if ((d2 > d1) && (d2 < rcut_oh)) {
                                adjmatrix[i][n_pro_n+n_ow+j] = 1;
                                conn_matrix[i][n_pro_n+n_ow+j] = 1;            //  i -> j is connected
								hb_h_ndx[i][n_pro_n+n_ow+j] = h_ndx[k]+1;
								// cout << network_ndx[i] << " " << network_ndx[n_pro_n+n_ow+j] <<  " " << hb_h_ndx[i][n_pro_n+n_ow+j] << " " << adjmatrix[i][n_pro_n+n_ow+j] << " " << adjmatrix[n_pro_n+n_ow+j][i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            //}else if((d2 < d1) && (d1 < rcut_oh) && (d2 < r_bond)) {
							}else if((d2 < d1) && (d1 < rcut_oh)){
                                adjmatrix[n_pro_n+n_ow+j][i] = 1;
                                conn_matrix[n_pro_n+n_ow+j][i] = 1;            //  j -> i is connected
								hb_h_ndx[n_pro_n+n_ow+j][i] = h_ndx[k]+1;
								// cout << network_ndx[i] << " " << network_ndx[n_pro_n+n_ow+j] <<  " " << hb_h_ndx[n_pro_n+n_ow+j][i] << " " << adjmatrix[i][n_pro_n+n_ow+j] << " " << adjmatrix[n_pro_n+n_ow+j][i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            //}else if((d2 == d1) && (d1 < r_bond)){
							}else if((d2 == d1) && (d1 < rcut_oh)){
								adjmatrix[i][n_pro_n+n_ow+j] = 1;
								adjmatrix[n_pro_n+n_ow+j][i] = 1;
                                conn_matrix[i][n_pro_n+n_ow+j] = 1;
                                conn_matrix[n_pro_n+n_ow+j][i] = 1;
								hb_h_ndx[i][n_pro_n+n_ow+j] = h_ndx[k]+1;
								hb_h_ndx[n_pro_n+n_ow+j][i] = h_ndx[k]+1;
								// cout << network_ndx[i] << " " << network_ndx[n_pro_n+n_ow+j] <<  " " << hb_h_ndx[i][n_pro_n+n_ow+j] << " " << adjmatrix[i][n_pro_n+n_ow+j] << " " << adjmatrix[n_pro_n+n_ow+j][i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            }
							
                        }
                        
                    }
                }
            }
		}
        for(int i=0; i<n_ow; ++i){
            for(int j=(i+1); j<n_ow; ++j){
                if(dist_calc(L, ow_pos[i], ow_pos[j]) < rcut_oo) { //distance between two oxygen atoms
                    for(int k=0; k<n_hw+n_pro_h; ++k){
                        hb_ang = angle_calc(L, ow_pos[i], h_pos[k], ow_pos[j]);
                        if (hb_ang > angcut_oho){
                            d1 = dist_calc(L, ow_pos[i], h_pos[k]);
                            d2 = dist_calc(L, ow_pos[j], h_pos[k]);
							//if ((d2 > d1) && (d2 < rcut_oh) && (d1 < r_bond)) {
                            if ((d2 > d1) && (d2 < rcut_oh)) {
                                adjmatrix[n_pro_n+i][n_pro_n+j] = 1;
                                conn_matrix[n_pro_n+i][n_pro_n+j] = 1;            //  i -> j is connected
								hb_h_ndx[n_pro_n+i][n_pro_n+j] = h_ndx[k]+1;
								// cout << network_ndx[n_pro_n+i] << " " << network_ndx[n_pro_n+j] <<  " " << hb_h_ndx[n_pro_n+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+j][n_pro_n+i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            //}else if((d2 < d1) && (d1 < rcut_oh) && (d2 < r_bond)) {
							}else if((d2 < d1) && (d1 < rcut_oh)){
                                adjmatrix[n_pro_n+j][n_pro_n+i] = 1;
                                conn_matrix[n_pro_n+j][n_pro_n+i] = 1;            //  j -> i is connected
								hb_h_ndx[n_pro_n+j][n_pro_n+i] = h_ndx[k]+1;
								// cout << network_ndx[n_pro_n+i] << " " << network_ndx[n_pro_n+j] <<  " " << hb_h_ndx[n_pro_n+j][n_pro_n+i] << " " << adjmatrix[n_pro_n+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+j][n_pro_n+i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            //}else if((d2 == d1) && (d1 < r_bond)){
							}else if((d2 == d1) && (d1 < rcut_oh)){
                                adjmatrix[n_pro_n+i][n_pro_n+j] = 1;
                                adjmatrix[n_pro_n+j][n_pro_n+i] = 1;
                                conn_matrix[n_pro_n+i][n_pro_n+j] = 1;            //  i -> j is connected
                                conn_matrix[n_pro_n+j][n_pro_n+i] = 1;            //  j -> i is connected
								hb_h_ndx[n_pro_n+i][n_pro_n+j] = h_ndx[k]+1;
								hb_h_ndx[n_pro_n+j][n_pro_n+i] = h_ndx[k]+1;
								// cout << network_ndx[n_pro_n+i] << " " << network_ndx[n_pro_n+j] <<  " " << hb_h_ndx[n_pro_n+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+j][n_pro_n+i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            }
							
                        }
                    }
                }
            }
        }
		for(int i = 0; i < n_pro_o; ++i){
			for(int j = 0; j < n_ow; ++j){
                if(dist_calc(L, pro_pos[pro_o_ndx[i]], ow_pos[j]) < rcut_oo) { //distance between two oxygen atoms
                    for(int k = 0; k < n_hw+n_pro_h; ++k){
                        hb_ang = angle_calc(L, pro_pos[pro_o_ndx[i]], h_pos[k], ow_pos[j]);
                        if (hb_ang > angcut_oho){
                            d1 = dist_calc(L, pro_pos[pro_o_ndx[i]], h_pos[k]);
                            d2 = dist_calc(L, ow_pos[j], h_pos[k]);
							//if ((d2 > d1) && (d2 < rcut_oh) && (d1 < r_bond)) {
                            if ((d2 > d1) && (d2 < rcut_oh)) {
                                adjmatrix[n_pro_n+n_ow+i][n_pro_n+j] = 1;
                                conn_matrix[n_pro_n+n_ow+i][n_pro_n+j] = 1;            //  i -> j is connected
								hb_h_ndx[n_pro_n+n_ow+i][n_pro_n+j] = h_ndx[k]+1;
								// cout << network_ndx[n_pro_n+n_ow+i] << " " << network_ndx[n_pro_n+j] <<  " " << hb_h_ndx[n_pro_n+n_ow+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+n_ow+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+j][n_pro_n+n_ow+i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            //}else if((d2 < d1) && (d1 < rcut_oh) && (d2 < r_bond)) {
							}else if((d2 < d1) && (d1 < rcut_oh)){
                                adjmatrix[n_pro_n+j][n_pro_n+n_ow+i] = 1;
                                conn_matrix[n_pro_n+j][n_pro_n+n_ow+i] = 1;            //  j -> i is connected
								hb_h_ndx[n_pro_n+j][n_pro_n+n_ow+i] = h_ndx[k]+1;
								// cout << network_ndx[n_pro_n+n_ow+i] << " " << network_ndx[n_pro_n+j] <<  " " << hb_h_ndx[n_pro_n+j] << " " << adjmatrix[n_pro_n+n_ow+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+j][n_pro_n+n_ow+i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            //}else if((d2 == d1) && (d1 < r_bond)){
							}else if((d2 == d1) && (d1 < rcut_oh)){
                                adjmatrix[n_pro_n+n_ow+i][n_pro_n+j] = 1;
                                adjmatrix[n_pro_n+j][n_pro_n+n_ow+i] = 1;
                                conn_matrix[n_pro_n+n_ow+i][n_pro_n+j] = 1;            //  i -> j is connected
                                conn_matrix[n_pro_n+j][n_pro_n+n_ow+i] = 1;            //  j -> i is connected
								hb_h_ndx[n_pro_n+n_ow+i][n_pro_n+j] = h_ndx[k]+1;
								hb_h_ndx[n_pro_n+j][n_pro_n+n_ow+i] = h_ndx[k]+1;
								// cout << network_ndx[n_pro_n+n_ow+i] << " " << network_ndx[n_pro_n+j] <<  " " << hb_h_ndx[n_pro_n+n_ow+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+n_ow+i][n_pro_n+j] << " " << adjmatrix[n_pro_n+j][n_pro_n+n_ow+i] << " " << d1 <<" " << d2 << " " << hb_ang << endl;
                            }
							
                        }
                        
                    }
                }
            }
		}
        /*~~~~~~~~~~~~~~~~~Dijkstra Alogorithm => Shortest Path ~~~~~~~~~~~~~~~~~~~~~~~*/
        length_file << n_frame << " ";
        path_file << n_frame << " ";
        in_file << n_frame << " ";
        out_file << n_frame << " ";
        
        int prev[V], visited[V], dist[V];
        int start, start1;
//        int length;
        int d, mc, u;
        int sum_in, sum_out;
        
        /*Initialization of path algorithm*/
        for(int j=0; j<V; j++){
            dist[j] = INF;
            visited[j] = 0;
            prev[j] =-1;
        }
                
        /*Dijkstra's algorithm*/
        start=s;
        dist[start] = 0;
        visited[start] = 1;
        do {
            mc=INF;     //!minimum cost between all paths
            u=1;
            for(int j=0; j<V; j++){
                if((j == start) || (visited[j] == 1)) continue;
                d = dist[start]+adjmatrix[start][j];
                if(d < dist[j]) {
                    dist[j] = d;
                    prev[j] = start;
                }
                if(mc > dist[j]) {
                    mc = dist[j];
                    u = j;
                }
            }
            if(mc == INF) break;
            start = u;
            visited[start] = 1;
        }while(visited[g]==0);
                
        /*Display*/
//        length=0;
        path_file << network_ndx[g] << " ";
        start1=g;
        int prev_start1;
        for(;;){
            prev_start1 = start1;
            start1=prev[start1];
//            cout << dist[g] << " " << adjmatrix[s][g] << endl;
//            cout << start1 << " " << ow_ndx[start1] << endl;
            if ((start1 == s) || (start1 == -1)) break;
//            length+=1;
            path_file << hb_h_ndx[start1][prev_start1] << " " << network_ndx[start1] << " ";
            sum_in=0;
            sum_out=0;
            for(int j=0; j<V; j++){
				sum_in += conn_matrix[j][start1];
				sum_out += conn_matrix[start1][j];
			}
            in_file << sum_in << " ";
            out_file << sum_out << " ";
        }
        if(start1==-1){
            path_file << start1 << endl;
            // proton_chk_file << "#" << endl;
        }
        else{
            path_file << hb_h_ndx[s][prev_start1] << " " << network_ndx[s] << endl;
            // proton_chk_file << conn_matrix[34][114] - conn_matrix[114][34] << " " << conn_matrix[53][34] - conn_matrix[34][53] << " " << conn_matrix[18][53] - conn_matrix[53][18] << " "<< endl;
            //oxygen index -> atom index: 114 -> 354, 34 -> 115, 53 -> 172, 18 -> 67
//            length+=1;
        }
        if (dist[g] == INF) dist[g]=0;
        length_file << dist[g] << endl;
        path_file << endl;
        in_file << endl;
        out_file << endl;
        n_frame +=1;
    }
    traj_file.close();
    sum_file.close();
    length_file.close();
    in_file.close();
    out_file.close();

    t_stop=clock();
    double runTime = (double)(t_stop-t_start) / (double)CLOCKS_PER_SEC;
    cout << "Total time: " << runTime << " [seconds] \n";

    delete[] ow_pos;
    delete[] h_pos;
    delete[] pro_pos;

    return 0;
}
/* distance between A--B */
double dist_calc(double L[][DIM], double *x,double *y){

    double temp;
    double z=0.;
    
    for(int i=0; i<DIM; ++i){
        temp = x[i]-y[i];         //x is one atom, y is another atom
        temp -= nearbyint(temp/L[i][i])*L[i][i];
        z += temp*temp;
    }
    z=sqrt(z);
    return z;
}

/* angle for D--H--A : angle between vector H-->D and vector H-->A */
double angle_calc(double L[][DIM], double *x, double *y, double *z){

    double t[3],u[3];
    double s,temp;
    
    for(int i=0; i<DIM; ++i){
        temp = x[i]-y[i];     //x is donor, y is hydrogen
        t[i] = temp - nearbyint(temp/L[i][i])*L[i][i];
        temp = z[i]-y[i];     //z is acceptor, y is hydrogen
        u[i] = temp - nearbyint(temp/L[i][i])*L[i][i];
    }

    s = scalar_product(t,u)/(norm(t)*norm(u));
    s = acos(s)*180/pi;
    
    return s;
}

/* scalar product between vector x and y */
double scalar_product(double *x, double *y){
    
    double s = 0.;

    for(int i=0; i<DIM; ++i){
        s += x[i]*y[i];
    }

    return s;
}

/* normalization of vector x */
double norm(double *x){
    
    double s=0.;

    for(int i=0; i<DIM; ++i){
        s += x[i]*x[i];
    }
    
    s=sqrt(s);

    return s;
}

void printUsage(string program)
{
    cout << endl;
    cout << "Usage:" << endl;
    cout << "  " << program << " traj.xyz out_file starting_node ending_node" << endl;
    cout << "Arguments:" << endl;
	cout << "   traj.xyz   (Input, input file for the analysis)" << endl;
	cout << "   starting_node    (Input, starting node for the proton wire, first N of NH: \"0\")" << endl;
    cout << "   ending_node      (Input, ending node of proton wire, second N of NH: \"1\")" << endl;
	cout << endl;
    return;
}
