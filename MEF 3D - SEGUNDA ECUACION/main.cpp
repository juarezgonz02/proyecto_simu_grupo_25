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

        /*
         * @example Correct usage mef.exe input_file [no file extension]
         */
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

        /*
            Read of .dat file saving data in Mesh M
        */
        read_input(filename, &M);

        M.report();

        /**
         *  @name Global / Acumulative values for FEM calculations
         */
        ///@{
        int num_nodes = M.get_quantity(NUM_NODES); //

        int num_elements = M.get_quantity(NUM_ELEMENTS);

        Matrix K(num_nodes, num_nodes), local_Ks[num_elements];

        Vector b(num_nodes), local_bs[num_elements];
        ///@}

        /**
         * Finite element method 
         * after solving result in a
         * Linear equation system with the form
         *
         * K*X = B
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
         * Then colects all local systems into a global system, 
         * it´s an assembly process
         * 
         * see mef_process.hpp -> create_local_systems() for more details
         */

        cout << "Creating local systems...\n\n";
        create_local_systems(local_Ks, local_bs, num_elements, &M);

        
        cout << "Performing Assembly...\n\n";
        /**
         * @brief Assembly all local_ks and local_bs into a GLOBAL K and GLOBAL B
         * 
         * - Local Matrices represent a element with 4 nodes
         * - Global Matrices represent the solution of the entire mesh 
         * 
         * - Each element has 4 nodes, but this nodes are shared between several elements, so, its necesary
         * assembly this elements, in order to relate the calculated values into a single node value in the 
         * global matrices  
         */
        assembly(&K, &b, local_Ks, local_bs, num_elements, &M);

        /**
         * @brief Apply boundary condition 
         * 
         * In MEF there are two types of conditions:
         * which are predefined values that serve as a starting point in the 
         * calculation of the elements and define the behavior at the boundaries
         * of the domain.         
         */
         
        cout << "Applying Neumann Boundary Conditions...\n\n";
        apply_neumann_boundary_conditions(&b, &M);

        cout << "Applying Dirichlet Boundary Conditions...\n\n";
        apply_dirichlet_boundary_conditions(&K, &b, &M);


        cout << "Solving global system...\n\n";
        Vector X(b.get_size()), X_full(num_nodes);
        solve_system(&K, &b, &X);

        cout << "Preparing results...\n\n";
        merge_results_with_dirichlet(&X, &X_full, num_nodes, &M);

        cout << "Writing output file...\n\n";
        write_output(filename, &T_full);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }

    return 0;
}
