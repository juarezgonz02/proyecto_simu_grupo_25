/**
 * @file mef_utilities/mef_process.hpp
 *
 * @author Enmanuel Amaya, MSc. 2023 *
 * @author Modified By: Grupo 09 - sec 01
 *
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
 *  [((1)/(3360*J))((B^T)(A^T)(A*B))][T1, T2, T3, T4]=(J/105)[1,1,1,1]
 *
 *  Where k is the thermal conductivity
 *  Where B =
 *      [-1, 1, 0, 0]
 *      [-1, 0, 1, 0]
 *      [-1, 0, 0, 1]
 *  Where J = is the jacobian calculated by the element´s nodes X, Y coordinates
 *  Where V = is the volume of the isoparametric piramid
 */
#include <cmath>
using namespace std;

float calculate_local_jacobian(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4)
{
    // 3D MEF CHANGE

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

void calculate_B(Matrix *B)
{
    // MEF 3D CHANGE #
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

void calculate_local_A(Matrix *A, float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4)
{
    // 3D MEF CHANGE

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

    K->set_size(4, 4);

    float x1 = M->get_element(element_id)->get_node1()->get_x_coordinate(), y1 = M->get_element(element_id)->get_node1()->get_y_coordinate(), z1 = M->get_element(element_id)->get_node1()->get_z_coordinate(),
          x2 = M->get_element(element_id)->get_node2()->get_x_coordinate(), y2 = M->get_element(element_id)->get_node2()->get_y_coordinate(), z2 = M->get_element(element_id)->get_node2()->get_z_coordinate(),
          x3 = M->get_element(element_id)->get_node3()->get_x_coordinate(), y3 = M->get_element(element_id)->get_node3()->get_y_coordinate(), z3 = M->get_element(element_id)->get_node3()->get_z_coordinate(),
          x4 = M->get_element(element_id)->get_node4()->get_x_coordinate(), y4 = M->get_element(element_id)->get_node4()->get_y_coordinate(), z4 = M->get_element(element_id)->get_node4()->get_z_coordinate();



    float J = calculate_local_jacobian(x1, y1, z1,
                                       x2, y2, z2,
                                       x3, y3, z3,
                                       x4, y4, z4);

    Matrix B(3, 4), A(3, 3);
    calculate_B(&B);
    calculate_local_A(&A, x1, y1, z1,
                      x2, y2, z2,
                      x3, y3, z3,
                      x4, y4, z4);
    // A.show();

    Matrix Bt(4, 3), At(3, 3);
    transpose(&B, 3, 4, &Bt);
    transpose(&A, 3, 3, &At);
     
     //Bt.show(); 
     
     //At.show();
    
    if(J == 0 || isnan(J)){
        J = 0.000006;
    }

    Matrix res1, res2, res3;
    
    product_matrix_by_matrix(&A, &B, &res1);

    product_matrix_by_matrix(&At, &res1, &res2);
    
    product_matrix_by_matrix(&Bt, &res2, &res3);
    
    product_scalar_by_matrix(1 / (3360 * J), &res3, 4, 4, K);
    
    // cout << "\t\tLocal matrix created for Element " << element_id+1 << ": "; K->show(); cout << "\n";

   // K->show();
}

void create_local_b(Vector *b,int element_id, Mesh *M)
{

    // MEF 3D CHANGE

    b->set_size(4);

    float x1 = M->get_element(element_id)->get_node1()->get_x_coordinate(), y1 = M->get_element(element_id)->get_node1()->get_y_coordinate(), z1 = M->get_element(element_id)->get_node1()->get_z_coordinate(),
          x2 = M->get_element(element_id)->get_node2()->get_x_coordinate(), y2 = M->get_element(element_id)->get_node2()->get_y_coordinate(), z2 = M->get_element(element_id)->get_node2()->get_z_coordinate(),
          x3 = M->get_element(element_id)->get_node3()->get_x_coordinate(), y3 = M->get_element(element_id)->get_node3()->get_y_coordinate(), z3 = M->get_element(element_id)->get_node3()->get_z_coordinate(),
          x4 = M->get_element(element_id)->get_node4()->get_x_coordinate(), y4 = M->get_element(element_id)->get_node4()->get_y_coordinate(), z4 = M->get_element(element_id)->get_node4()->get_z_coordinate();
    float J = calculate_local_jacobian(x1, y1, z1,
                                       x2, y2, z2,
                                       x3, y3, z3,
                                       x4, y4, z4);

    b->set(J / 105, 0);
    b->set(J / 105, 1);
    b->set(J / 105, 2);
    b->set(J / 105, 3);

    // cout << "\t\tLocal vector created for Element " << element_id+1 << ": "; b->show(); cout << "\n";
}

void create_local_systems(Matrix *Ks, Vector *bs,int num_elements, Mesh *M)
{

    for (int e = 0; e < num_elements; e++)
    {

        cout << "\tCreating local system for Element " << e + 1 << "...\n\n";

        create_local_K(&Ks[e], e, M);
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
    K->init();
    b->init();
    // K->show(); b->show();

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
        // cout << "\t\t"; K->show(); cout << "\t\t"; b->show(); cout << "\n";
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

void merge_results_with_dirichlet(Vector *X, Vector *Tf, int n, Mesh *M)
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
            Tf->set(X->get(cont_T), i);
            cont_T++;
        }
    }
}

void solve_system(Matrix *K, Vector *b, Vector *X)
{
    int n = K->get_nrows();

    Matrix Kinv(n, n);
    cout << n << "\n\n\n";
    //K->show();
    cout << "\tCalculating inverse of global matrix K...\n\n";
    calculate_inverse(K, n, &Kinv);

    //Kinv.show();
    cout << "\tPerforming final calculation...\n\n";
    product_matrix_by_vector(&Kinv, b, n, n, X);
}
