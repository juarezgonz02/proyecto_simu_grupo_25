/**
 * @file input_output.hpp
 * @author Enmanuel Amaya, MSc. 2023
 * @author Modified By: Grupo 09 - sec 01
 * @brief File manipulation
 * 
 * Input file must be generated with this format
 * 
    [k] [Q]
    [Dirichlet condtion value] [Neuman condition value]
    [num_nodes] [num_elements] [num_dirichlet_nodes] [num_nueman_nodes]

    Coordinates
    ...
    [node_id][x][y][z]
    ...
    EndCoordinates

    Elements
    ...
    [element_id][node_1_id][node_2_id][node_3_id][node_4_id]
    ...
    EndElements

    Dirichlet
    ...
    [nodes_with_dirichlet_id]
    ...
    EndDirichlet

    Neumann
    ...
    [nodes_with_nuemann_id]
    ...
    EndNeumann 
 */

#include <fstream>
#include <iostream>
#include <string>
using namespace std;

void read_input(string filename, Mesh* M){
    string line;
    float k, Q, T_bar, T_hat;
    int num_nodes, num_elements, num_dirichlet, num_neumann;
    ifstream dat_file(filename+".dat");

    dat_file >> k >> Q >> T_bar >> T_hat >> num_nodes >> num_elements >> num_dirichlet >> num_neumann;
    
    M->set_problem_data(k,Q);
    M->set_quantities(num_nodes,num_elements,num_dirichlet,num_neumann);

    M->init_arrays();

    // While reading there are lines that limit sections, such as Elements / Conditions / Nodes  
    // So read that line, saving a temp variable, in order to skip the lines 
    dat_file >> line;


    /**
     * @brief Read node list
     * 
     */
    for(int i = 0; i < num_nodes; i++){
        int id;
        float x, y, z;
        dat_file >> id >> x >> y >> z;

        M->insert_node(new Node(id,x,y, z), i);

    }

     /**
     * @brief Read element list
     * 
     * Skip section limit lines and read element´s nodes and id
     */
    dat_file >> line >> line;
    cout << "num_elements " << num_elements<< "\n";
    for(int i = 0; i < num_elements; i++){

        /**
         * @name 3D MEF CHANGE 
         * 
         * Reading and Adding 4th node from input file to Mesh  
         */
        int id, node1_id, node2_id, node3_id, node4_id;

        dat_file >> id >> node1_id >> node2_id >> node3_id >> node4_id ;
        
        M->insert_element(new Element(id, M->get_node(node1_id-1), M->get_node(node2_id-1), M->get_node(node3_id-1), M->get_node(node4_id-1)), i);
    }


     /**
     * @brief Read conditions list
     * 
     * Skip section limit lines and read nodes with dirichlet conditions 
     */
    dat_file >> line >> line;

    for(int i = 0; i < num_dirichlet; i++){
       int id;
        dat_file >> id;
        M->insert_dirichlet_condition(new Condition(M->get_node(id-1), T_bar), i);
    }

    dat_file >> line >> line;

     /**
     * @brief Read conditions list
     * 
     * Skip section limit lines and read nodes with neumman conditions 
     */
    for(int i = 0; i < num_neumann; i++){
       int id;
        dat_file >> id;

        M->insert_neumann_condition(new Condition(M->get_node(id-1), T_hat), i);
    }

    //ALWAYS CLOSE READ FILE STREAM 
    dat_file.close();
}

/**
 * Output file must have the following format
 *
    GiD Post Results File 1.0
    Result "Temperature" "Load Case 1" 1 Scalar OnNodes
    ComponentNames "T"
    Values
    ...
    [node_id]     [node_value]
    ...
    End values
 * 
 */

/**
 * @brief Output Writter
 *  
 */
void write_output(string filename, Vector* T){

    /**
     * 
     */
    ofstream res_file(filename+".post.res");
    
    res_file << "GiD Post Results File 1.0\n";

    int n = T->get_size();

    res_file << "Result \"Temperature\" \"Load Case 1\" " << 1 << " Scalar OnNodes\n";
    res_file << "ComponentNames \"T\"\n";
    res_file << "Values\n";

    for(int i = 0; i < n; i++)
        res_file << i+1 << "     " << T->get(i) << "\n";

    res_file << "End values\n";

    res_file.close();
}
