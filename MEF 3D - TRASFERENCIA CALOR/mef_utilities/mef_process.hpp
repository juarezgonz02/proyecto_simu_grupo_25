/**
 * @file mef_utilities/mef_process.hpp
 *
 * @author Enmanuel Amaya, MSc. 2023 *
 * @author Modified By: Grupo 09 - sec 01
 *
 * @brief 
 * @version 1
 * @date 2023-06-30
 *
 * C++ FEM process utils
 *
 *  After having performed the development of integrals and interpolation
 *  as a result of the FEM we obtain a series of general equations and matrices
 *  for any number of nodes, elements, nodes with boundary conditions.
 *  This document contains the implementation of those equations and matrix
 *  manipulation in c++.
 *
 *  Especial implementation of FEM 3D for Heat Transfer model has the form of
 *
 *  [((k*V)/(J*J))((B^T)(A^T)(A*B))][T1, T2, T3, T4]=(Q*J/24)[1,1,1,1]
 *
 *  Where k is the thermal conductivity
 *  Where B =
 *      [-1, 1, 0, 0]
 *      [-1, 0, 1, 0]
 *      [-1, 0, 0, 1]
 *  Where J = is the jacobian calculated by the element´s nodes X, Y, Z coordinates
 *  Where V = is the volume of the isoparametric piramid
 */
#include <cmath>
using namespace std;

float calculate_local_volume(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4);
float calculate_local_volume(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4)
{
    // 3D MEF CHANGE
    Matrix volume_matrix(3  , 3);
    volume_matrix.set((x2 - x1), 0, 0);
    volume_matrix.set((y2 - y1), 0, 1);
    volume_matrix.set((z2 - z1), 0, 2);

    volume_matrix.set((x3 - x1), 1, 0);
    volume_matrix.set((y3 - y1), 1, 1);
    volume_matrix.set((z3 - z1), 1, 2);

    volume_matrix.set((x4 - x1), 2, 0);
    volume_matrix.set((y4 - y1), 2, 1);
    volume_matrix.set((z4 - z1), 2, 2);

    return (1.0 / 6.0) * abs(determinant(&volume_matrix));
}

/**
 * @brief Calculates the jacobian for a element based in
 * its coordinates
 * 
 * The Jacobian in the FEM is a matrix used to map the physical coordinates of a finite element
 * to its natural coordinates, which allows the formulation and solution of the equations of the 
 * problem in the reference domain.
 * 
 * @param x1 Coordinate X for node 1    
 * @param y1 Coordinate Y for node 1  
 * @param z1 Coordinate Z for node 1  
 * @param x2 Coordinate X for node 2    
 * @param y2 Coordinate Y for node 2  
 * @param z2 Coordinate Z for node 2  
 * @param x3 Coordinate X for node 3    
 * @param y3 Coordinate Y for node 3  
 * @param z3 Coordinate Z for node 3  
 * @param x4 Coordinate X for node 4    
 * @param y4 Coordinate Y for node 4  
 * @param z4 Coordinate Z for node 4  
 * @return float 
 */

float calculate_local_jacobian(float x1 , float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4)
{
    // 3D MEF CHANGE
    /**
     * Jacobian matrix has the form of 
     *     [X2 - X1, X3 - X1, X4 - X1]
     * J = [Y2 - Y1, Y3 - Y1, Y4 - Y1]
     *     [Z2 - Z1, Z3 - Z1, Z4 - Z1]
     * 
     * Then what we need is the determinant of this matrix
     */
    Matrix jacobian_matrix(3, 3);   
    jacobian_matrix.set((x2 - x1), 0, 0);
    jacobian_matrix.set((x3 - x1), 0, 1);
    jacobian_matrix.set((x4 - x1), 0, 2);

    jacobian_matrix.set((y2 - y1), 1, 0);
    jacobian_matrix.set((y3 - y1), 1, 1);
    jacobian_matrix.set((y4 - y1), 1, 2);

    jacobian_matrix.set((z2 - z1), 2, 0);
    jacobian_matrix.set((z3 - z1), 2, 1);
    jacobian_matrix.set((z4 - z1), 2, 2);

    return determinant(&jacobian_matrix);
}

/**
 * @brief Fills a matrix with predefined values
 * 
 * @param B Reference to the matrix to be filled, has to have the form of 4x4
 */
void calculate_B(Matrix *B)
{
    // MEF 3D CHANGE #

    /*
    * Matrix has the form of 
    * B =
    *    [-1, 1, 0, 0]
    *    [-1, 0, 1, 0]
    *    [-1, 0, 0, 1]
    */

    B->set(-1, 0, 0);
    B->set(1, 0, 1);
    B->set(0, 0, 2);
    B->set(0, 0, 3);

    B->set(-1, 1, 0);
    B->set(0, 1, 1);
    B->set(1, 1, 2);
    B->set(0, 1, 3);

    B->set(-1, 2, 0);
    B->set(0, 2, 1);
    B->set(0, 2, 2);
    B->set(1, 2, 3);
}

/**
 * @brief Calculates the *attached* Matrix to the Jacobian for a element based in
 * its coordinates
 * 
 * The Jacobian in the FEM is a matrix used to map the physical coordinates of a finite element
 * to its natural coordinates, which allows the formulation and solution of the equations of the 
 * problem in the reference domain.
 * 
 * @param A Takes a matrix by reference to fill with the values of the attached
 * @param x1 Coordinate X for node 1    
 * @param y1 Coordinate Y for node 1  
 * @param z1 Coordinate Z for node 1  
 * @param x2 Coordinate X for node 2    
 * @param y2 Coordinate Y for node 2  
 * @param z2 Coordinate Z for node 2  
 * @param x3 Coordinate X for node 3    
 * @param y3 Coordinate Y for node 3  
 * @param z3 Coordinate Z for node 3  
 * @param x4 Coordinate X for node 4    
 * @param y4 Coordinate Y for node 4  
 * @param z4 Coordinate Z for node 4  
 * @return float 
 */

void calculate_local_A(Matrix *A, float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4)
{


    A->set((y3 - y1) * (z4 - z1) - (y4 - y1) * (z3 - z1), 0, 0);
    A->set(-(x3 - x1) * (z4 - z1) + (x4 - x1) * (z3 - z1), 0, 1);
    A->set((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1), 0, 2);

    A->set(-(y2 - y1) * (z4 - z1) + (y4 - y1) * (z2 - z1), 1, 0);
    A->set((x2 - x1) * (y4 - y1) + (x4 - x1) * (y2 - y1), 1, 1);
    A->set(-(x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1), 1, 2);

    A->set((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1), 2, 0);
    A->set(-(x2 - x1) * (z3 - z1) + (x3 - x1) * (z2 - z1), 2, 1);
    A->set((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1), 2, 2);
}

void create_local_K(Matrix *K,int element_id, Mesh *M)
{
    // MEF 3D CHANGE
    /**
     * Initializate each local_k matrix, 
     * Setted size at 4x4, since a element in this MEF aplication is a tetahedron, so it uses 4 nodes
     */
    K->set_size(4, 4);


    // Get element´s group of coordinates
    float k = M->get_problem_data(THERMAL_CONDUCTIVITY);
    float x1 = M->get_element(element_id)->get_node1()->get_x_coordinate(), y1 = M->get_element(element_id)->get_node1()->get_y_coordinate(), z1 = M->get_element(element_id)->get_node1()->get_z_coordinate(),
          x2 = M->get_element(element_id)->get_node2()->get_x_coordinate(), y2 = M->get_element(element_id)->get_node2()->get_y_coordinate(), z2 = M->get_element(element_id)->get_node2()->get_z_coordinate(),
          x3 = M->get_element(element_id)->get_node3()->get_x_coordinate(), y3 = M->get_element(element_id)->get_node3()->get_y_coordinate(), z3 = M->get_element(element_id)->get_node3()->get_z_coordinate(),
          x4 = M->get_element(element_id)->get_node4()->get_x_coordinate(), y4 = M->get_element(element_id)->get_node4()->get_y_coordinate(), z4 = M->get_element(element_id)->get_node4()->get_z_coordinate();


    // Calculate element volumen using local coordinates
    float volume = calculate_local_volume(x1, y1, z1,
                                          x2, y2, z2,
                                          x3, y3, z3,
                                          x4, y4, z4);


    // Calculate element jacobian using local coordinates
    float J = calculate_local_jacobian(x1, y1, z1,
                                       x2, y2, z2,
                                       x3, y3, z3,
                                       x4, y4, z4);

    /**
     * @brief Initializate B and A matrices  
     * 
     * These matrices are part of the result of solving the integrals of the FEM process.
     * 
     * *For efficiency reasons, all matrices are always passed by reference.
     * 
     * Matrix B has 3 column, 4 rows form
     * Matrix A is a 3x3 matrix
     */
    Matrix B(3, 4), A(3, 3);
    
    calculate_B(&B);
    calculate_local_A(&A, x1, y1, z1,
                      x2, y2, z2,
                      x3, y3, z3,
                      x4, y4, z4);
    // A.show();

    /**
     * @brief Creates a transposed matrix by exchanging the dimensions of the original matrix.  
     *  
     */

    Matrix Bt(4, 3), At(3, 3);
    transpose(&B, 3, 4, &Bt);
    transpose(&A, 3, 3, &At);
         
    /**
     * @brief Dump Patch 
     * There´s some combination of coordinates that generates 
     * a jacobian or a volume with value 0 
     */
    if(J == 0 || isnan(J)){
        J = 0.000006;
    }
    if(volume == 0 || isnan(volume)){
        volume = 0.000006;
    }

    /**
     * @brief Multiply matrix
     * 
     * Remembering the form of the local K, we must perform the multiplications
     * from right to left, since the multiplication of matrices is not commutative,
     * but associative, so we take two matrices and the result is multiplied by the following one

     * [((k*V)/(J*J))((B^T)(A^T)(A*B))]  
     * 
     */
    Matrix res1, res2, res3;

    //Multiply A X B
    product_matrix_by_matrix(&A, &B, &res1);

    //Multiply (A^T)(A X B)
    product_matrix_by_matrix(&At, &res1, &res2);
    
    //Multiply (B^T)[(A^T)(A X B)]
    product_matrix_by_matrix(&Bt, &res2, &res3);
    
    //Multiply escalar value (k*v/jj ) and result matrix [(B^T)[(A^T)(A X B)]]
    product_scalar_by_matrix(k * volume / (J * J), &res3, 4, 4, K);
    
}

void create_local_b(Vector *b,int element_id, Mesh *M)
{

    // MEF 3D CHANGE
    /**
     * Initializate each local_b vector, 
     * Setted size at 4, since a element in this MEF aplication is a tetahedron, so it uses 4 nodes
     */
    b->set_size(4);

    //Problem data    
    float Q = M->get_problem_data(HEAT_SOURCE);

    // Get element´s group of coordinates
    float x1 = M->get_element(element_id)->get_node1()->get_x_coordinate(), y1 = M->get_element(element_id)->get_node1()->get_y_coordinate(), z1 = M->get_element(element_id)->get_node1()->get_z_coordinate(),
          x2 = M->get_element(element_id)->get_node2()->get_x_coordinate(), y2 = M->get_element(element_id)->get_node2()->get_y_coordinate(), z2 = M->get_element(element_id)->get_node2()->get_z_coordinate(),
          x3 = M->get_element(element_id)->get_node3()->get_x_coordinate(), y3 = M->get_element(element_id)->get_node3()->get_y_coordinate(), z3 = M->get_element(element_id)->get_node3()->get_z_coordinate(),
          x4 = M->get_element(element_id)->get_node4()->get_x_coordinate(), y4 = M->get_element(element_id)->get_node4()->get_y_coordinate(), z4 = M->get_element(element_id)->get_node4()->get_z_coordinate();
    

    // Calculate jacobian using local coordinates
    float J = calculate_local_jacobian(x1, y1, z1,
                                       x2, y2, z2,
                                       x3, y3, z3,
                                       x4, y4, z4);

    // local_b has the form of 
    /*
        (QJ/24)[1]
        (QJ/24)[1]
        (QJ/24)[1]
        (QJ/24)[1]
    */
    b->set(Q * J / 24, 0);
    b->set(Q * J / 24, 1);
    b->set(Q * J / 24, 2);
    b->set(Q * J / 24, 3);

}

void create_local_systems(Matrix *Ks, Vector *bs,int num_elements, Mesh *M)
{

    //Creates the local system for each element 
    for (int e = 0; e < num_elements; e++)
    {

        cout << "\tCreating local system for Element " << e + 1 << "...\n\n";
        /**
         * @brief Create a local K object
         * Creating the local K means calculate each piece of the formula for this element
         * [((k*V)/(J*J))((B^T)(A^T)(A*B))] 
         * and multiply
         * the resulting matrices
         */
        create_local_K(&Ks[e], e, M);
        
        /**
         * @brief Create a local b object
         * Creating the local b means calculate the jacobian (J) for each element and replace 
         * it in the formula (Q*J/24)[1,1,1,1]
         */
        create_local_b(&bs[e], e, M);
    }
}


void assembly_K(Matrix *K, Matrix *local_K,int index1,int index2,int index3,int index4)
{

    // 3D MEF CHANGE
    K->add(local_K->get(0, 0), index1, index1);
    K->add(local_K->get(0, 1), index1, index2);
    K->add(local_K->get(0, 2), index1, index3);
    K->add(local_K->get(0, 3), index1, index4);

    K->add(local_K->get(1, 0), index2, index1);
    K->add(local_K->get(1, 1), index2, index2);
    K->add(local_K->get(1, 2), index2, index3);
    K->add(local_K->get(1, 3), index2, index4);

    K->add(local_K->get(2, 0), index3, index1);
    K->add(local_K->get(2, 1), index3, index2);
    K->add(local_K->get(2, 2), index3, index3);
    K->add(local_K->get(2, 3), index3, index4);

    K->add(local_K->get(3, 0), index4, index1);
    K->add(local_K->get(3, 1), index4, index2);
    K->add(local_K->get(3, 2), index4, index3);
    K->add(local_K->get(3, 3), index4, index4);

    // local_K->show();
}

void assembly_b(Vector *b, Vector *local_b,int index1,int index2,int index3,int index4)
{

    // 3D MEF CHANGE
    b->add(local_b->get(0), index1);
    b->add(local_b->get(1), index2);
    b->add(local_b->get(2), index3);
    b->add(local_b->get(3), index4);
}

void assembly(Matrix *K, Vector *b, Matrix *Ks, Vector *bs,int num_elements, Mesh *M)
{

    /**
     * @brief Initializate values of matrix / vector filling with zeros 
     * 
     */
    K->init();
    b->init();


    //For each element
    for (int e = 0; e < num_elements; e++)
    {
        // 3D MEF CHANGE
        cout << "\tAssembling for Element " << e + 1 << "...\n\n";
       int index1 = M->get_element(e)->get_node1()->get_ID() - 1;
       int index2 = M->get_element(e)->get_node2()->get_ID() - 1;
       int index3 = M->get_element(e)->get_node3()->get_ID() - 1;
       int index4 = M->get_element(e)->get_node4()->get_ID() - 1;

        assembly_K(K, &Ks[e], index1, index2, index3, index4);
        assembly_b(b, &bs[e], index1, index2, index3, index4);
    }
}

void apply_neumann_boundary_conditions(Vector *b, Mesh *M)
{
   int num_conditions = M->get_quantity(NUM_NEUMANN);

    for (int c = 0; c < num_conditions; c++)
    {
        Condition *cond = M->get_neumann_condition(c);

       int index = cond->get_node()->get_ID() - 1;
        b->add(cond->get_value(), index);
    }
    /// cout << "\t\t"; b->show(); cout << "\n";
}

void add_column_to_RHS(Matrix *K, Vector *b, int col, float T_bar)
{
    for (int r = 0; r < K->get_nrows(); r++)
        b->add(-T_bar * K->get(r, col), r);
}

void apply_dirichlet_boundary_conditions(Matrix *K, Vector *b, Mesh *M)
{
   int num_conditions = M->get_quantity(NUM_DIRICHLET);
    int previous_removed = 0;

    for (int c = 0; c < num_conditions; c++)
    {
        Condition *cond = M->get_dirichlet_condition(c);

       int index = cond->get_node()->get_ID() - 1 - previous_removed;
        float cond_value = cond->get_value();

        K->remove_row(index);
        // K->show();
        //  b->show();
        b->remove_row(index);
        // b->show();

        add_column_to_RHS(K, b, index, cond_value);
        // b->show();

        K->remove_column(index);
        // K->show();

        previous_removed++;
    }
}

void merge_results_with_dirichlet(Vector *T, Vector *Tf, int n, Mesh *M)
{
    int num_dirichlet = M->get_quantity(NUM_DIRICHLET);

    int cont_dirichlet = 0;
    int cont_T = 0;
    for (int i = 0; i < n; i++)
    {
        if (M->does_node_have_dirichlet_condition(i + 1))
        {
            Condition *cond = M->get_dirichlet_condition(cont_dirichlet);
            cont_dirichlet++;

            float cond_value = cond->get_value();

            Tf->set(cond_value, i);
        }
        else
        {
            Tf->set(T->get(cont_T), i);
            cont_T++;
        }
    }
}

void solve_system(Matrix *K, Vector *b, Vector *T)
{
    int n = K->get_nrows();

    Matrix Kinv(n, n);
    cout << n << "\n\n\n";
    //K->show();
    cout << "\tCalculating inverse of global matrix K...\n\n";
    calculate_inverse(K, n, &Kinv);

    //Kinv.show();
    cout << "\tPerforming final calculation...\n\n";
    product_matrix_by_vector(&Kinv, b, n, n, T);
}
