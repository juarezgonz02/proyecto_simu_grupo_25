#include <iostream>
#include <cstdlib>

using namespace std;

#include "geometry/mesh.hpp"
#include "math_utilities/matrix_operations.hpp"
#include "mef_utilities/mef_process.hpp"
#include "gid/input_output.hpp"
/*
 * @brief MEF 3D
 *
 * Implementation for the Finite Element Method for a 3D mesh using data
 * from GID Mesh problemType generation
 *
 * @param argc number of params passed at execution
 * @param argv pointer to char chain witch represents every value passed at execution
 *
 */
int main(int argc, char **argv)
{

    try
    {
    
    
    // @example Correct usage mef.exe input_file [no file extension] 
    if (argc != 2)
    {
        cout << "Incorrect use of the program, it must be: mef filename\n";
        exit(EXIT_FAILURE);
    }
    

    /*
    Mesh representation declarations 
    */
    Mesh M;

    cout << "Reading geometry and mesh data...\n\n";

    /*
     Using string constructor from char* to string
    */
    string filename(argv[1]);
    //string filename = "input";
    /*
        Read of .dat file saving data in Mesh M
    */
    read_input(filename, &M);

    M.report();


    /**
     *  @name Global / Acumulative values for FEM calculations 
     *
     *  They are required by default to cooperate with the global C
     *  library's `FILE` streams, and to be available during program
     *  startup and termination. For more information, see the section of the
     *  manual linked to above.
     */
    ///@{
    int num_nodes = M.get_quantity(NUM_NODES); // 

    int num_elements = M.get_quantity(NUM_ELEMENTS);

    Matrix K(num_nodes, num_nodes), local_Ks[num_elements];

    Vector b(num_nodes), local_bs[num_elements];
    ///@}

    /**
     * Finite element method result in a 
     * Linear equation system with the form
     * 
     * KX = B
     * 
     * [a... b][x1]=[e]
     * [.... .][x2]=[f]
     * [c... d][x3]=[g]
     * 
     * Where K is N*N Matrix continaing the coefficients of the unknowns of the system
     * Where X is N*1 VECTOR continaing the unknowns of the system
     * 
     * Where B is N*1 VECTOR containing the result of each equation
     * 
     * This section calculates this Equation System for each element in the mesh
     * its called a local system. each local system is saved into local_Bs and local_Ks vectors 
     * 
     * Then colects all local systems into a global system, it´s an assembly process 
     */

    cout << "Creating local systems...\n\n";
    create_local_systems(local_Ks, local_bs, num_elements, &M);

    cout << "Performing Assembly...\n\n";
    assembly(&K, &b, local_Ks, local_bs, num_elements, &M);

    //K.show(); 
    //b.show();

    cout << "Applying Neumann Boundary Conditions...\n\n";
    apply_neumann_boundary_conditions(&b, &M);
    
    cout << "Applying Dirichlet Boundary Conditions...\n\n";
    apply_dirichlet_boundary_conditions(&K, &b, &M); 

    cout << "Solving global system...\n\n";
    Vector T(b.get_size()), T_full(num_nodes);
    solve_system(&K, &b, &T);

    cout << "Preparing results...\n\n";
    merge_results_with_dirichlet(&T, &T_full, num_nodes, &M);

    cout << "Writing output file...\n\n";
    write_output(filename, &T_full);
    
}
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    return 0;
}
