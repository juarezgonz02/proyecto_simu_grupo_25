/**
 * @file math_utilities/matrix_operations.hpp
 *
 * @author Enmanuel Amaya, MSc. 2023 *
 * @author Modified By: Grupo 09 - sec 01
 *
 * @brief Matrix Operations
 * @version 1
 * @date 2023-06-30
 *
 */

#include <cmath>
#include "vector.hpp"
#include "matrix.hpp"

/**
 * @brief Calculates the product of a matrix and a scalar
 *
 * A scalar value is a single value with no dimension, so takes each value in matrix,
 * and multiply by the scalar value
 *
 * @param n Rows of Matrix M
 * @param m Columns of matrix M
 */
void product_scalar_by_matrix(float scalar, Matrix *M, int n, int m, Matrix *R)
{
    for (int r = 0; r < n; r++)
        for (int c = 0; c < m; c++)
            R->set(scalar * M->get(r, c), r, c);
}

/**
 * @brief Performs matrix-vector multiplication.
 *
 * @param n Rows of Matrix M
 * @param m Columns of matrix M
 * @param R Output Vector
 */
void product_matrix_by_vector(Matrix *M, Vector *V, int n, int m, Vector *R)
{
    for (int r = 0; r < n; r++)
    {
        float acc = 0;
        for (int c = 0; c < n; c++)
            acc += M->get(r, c) * V->get(c);
        R->set(acc, r);
    }
}

/**
 * @brief Performs the multiplication of two matrices.
 *
 * It retrieves the dimensions of matrices A and B.
 * It checks if the number of columns in A is equal to the number of rows in B.
 * If not, it prints an error message and exits the program.
 *
 * If the dimensions are compatible, it initializes matrix R with the appropriate
 * size and initializes its values filling with zeros.
 *
 * It iterates through each row of A and each column of B to calculate the elements of R.
 *
 * For each element in R, it performs a nested loop that multiplies
 * corresponding elements from A and B, accumulating the sum of these products.
 *
 * The result of each sum is stored in the corresponding position of R.
 *
 */

void product_matrix_by_matrix(Matrix *A, Matrix *B, Matrix *R)
{
    int n = A->get_nrows(),
        m = A->get_ncols(),

        p = B->get_nrows(),
        q = B->get_ncols();

    if (m == p)
    {
        R->set_size(n, q);
        R->init();

        /**
         * loop throght each value in Matrix R
         */
        for (int r = 0; r < n; r++)
            for (int c = 0; c < q; c++)
                for (int i = 0; i < m; i++)
                    /**
                        - A->get(r, i) retrieves the element at the r-th row and i-th column of matrix A.
                          B->get(i, c) retrieves the element at the i-th row and c-th column of matrix B.

                        - Multply this values
                        - The resulting product is then added to the element at position (r, c)
                          in matrix R
                     */
                    R->add(A->get(r, i) * B->get(i, c), r, c);
    }
    else
    {
        cout << "Incompatibilidad de dimensiones al multiplicar matrices.\n\nAbortando...\n";
        exit(EXIT_FAILURE);
    }
}

float determinant(Matrix *M);

float determinant_auxiliar(Matrix *M)
{
    int n = M->get_ncols();
    float acc = 0;

    for (int c = 0; c < n; c++)
    {
        Matrix clon(n, n);
        M->clone(&clon);
        clon.remove_row(0);
        clon.remove_column(c);
        acc += pow(-1, c) * M->get(0, c) * determinant(&clon);
    }

    return acc;
}

float determinant(Matrix *M)
{
    float ans;
    switch (M->get_ncols())
    {
    case 1:
        ans = M->get(0, 0);
        break;
    case 2:
        ans = M->get(0, 0) * M->get(1, 1) - M->get(0, 1) * M->get(1, 0);
        break;
    case 3:
        ans = M->get(0, 0) * M->get(1, 1) * M->get(2, 2) - M->get(0, 0) * M->get(1, 2) * M->get(2, 1) - M->get(0, 1) * M->get(1, 0) * M->get(2, 2) + M->get(0, 1) * M->get(1, 2) * M->get(2, 0) + M->get(0, 2) * M->get(1, 0) * M->get(2, 1) - M->get(0, 2) * M->get(1, 1) * M->get(2, 0);
        break;
    default:
        ans = determinant_auxiliar(M);
    }
    return ans;
}

float get_minor(Matrix *M, int n, int r, int c)
{
    Matrix clon(n, n);
    M->clone(&clon);

    clon.remove_row(r);
    clon.remove_column(c);

    return determinant(&clon);
}

void conjugate_matrix(Matrix *M, int n, Matrix *C)
{
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++)
            C->set(pow(-1, r + c) * get_minor(M, n, r, c), r, c);
}

/**
 * @brief Calculates the transpose of a matrix
 *
 * The transpose of a matrix is an operation involving the exchange of rows
 * for columns in a given matrix. Given a matrix A of size m x n, the transpose
 * of A, denoted A^T, is a new matrix of size n x m where the elements of the rows
 * in A become the corresponding columns in A^T.
 *
 * Each value in the original matrix with position (row, column) will be set to the
 * transponse matrix to position (column, row)
 * @param M
 * @param n
 * @param m
 * @param T
 */
void transpose(Matrix *M, int n, int m, Matrix *T)
{
    for (int r = 0; r < n; r++)
        for (int c = 0; c < m; c++)
            T->set(M->get(r, c), c, r);
}

/**
 * @brief Implementation of CHOLESKY METHOD for inverse matrix computation 
 * @author Enmanuel Amaya, MSc. 2023 *
 *
 * See more at http://funes.uniandes.edu.co/8037/1/Alpizar2013Factorizacion.pdf
 * @param A Input matrix
 * @param n Matrix size
 * @param X Output Matrix
 */
void calculate_inverse(Matrix *A, int n, Matrix *X)
{
    Matrix L(n, n), Y(n, n);
    float acum;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                acum = 0;
                for (int k = 0; k < j; k++)
                {
                    acum += pow(L.get(j, k), 2);
                }

                if ((A->get(j, j) - acum) <= 0)
                {
                    L.set(0.000006, j, j);
                }
                else
                {
                    L.set(sqrt(A->get(j, j) - acum), j, j);
                }
            }
            else
            {
                if (i > j)
                {
                    acum = 0;
                    for (int k = 0; k < j; k++)
                    {
                        acum += L.get(i, k) * L.get(j, k);
                    }

                    L.set((1 / L.get(j, j)) * (A->get(i, j) - acum), i, j);
                }
                else
                {
                    L.set(0, i, j);
                }
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                Y.set(1 / L.get(i, i), i, i);
            }
            else
            {
                if (i > j)
                {
                    acum = 0;
                    for (int k = j; k < i; k++)
                    {
                        acum += L.get(i, k) * Y.get(k, j);
                    }
                    Y.set(-(1 / L.get(i, i)) * acum, i, j);
                }
                else
                {
                    Y.set(0, i, j);
                }
            }
        }
    }

    for (int i = n - 1; i >= 0; i--)
    {
        for (int j = 0; j < n; j++)
        {

            acum = 0;

            for (int k = i + 1; k < n; k++)
            {
                acum += L.get(k, i) * X->get(k, j);
            }

            X->set((1 / L.get(i, i)) * (Y.get(i, j) - acum), i, j);
        }
    }
}
