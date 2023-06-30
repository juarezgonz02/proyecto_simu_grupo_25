#include "heap.hpp"
#include "element.hpp"
#include "condition.hpp"

/**
 * @brief Heat Transfer Model Constants
 *
 * Declared for fast access throught problem_data array list
 */
enum parameter
{
    THERMAL_CONDUCTIVITY, // K in Heat Transfer Ecuation
    HEAT_SOURCE           // Q in Heat Transfer Ecuation
};

/**
 * @brief FEM Values sizes
 *
 * Declared for fast access throught quantities array list
 */
enum quantity
{
    NUM_NODES,     // nodes list size position
    NUM_ELEMENTS,  // elements list size position
    NUM_DIRICHLET, // dirichlet_conditions list size position
    NUM_NEUMANN    // neuman_conditions list size position
};

/**
 * @file geometry/mesh
 * @brief 3D MESH Modeling for FEM Calculations
 *
 * This implementation has as problem
 */
class Mesh
{

private:
    float problem_data[2]; // Array list of constants values for especified model
    int quantities[4];   // Array list of FEM values sizes

    /**
     *  @name FEM VALUES Colections
     *
     *  This values in this secction are a list of pointer colecctions
     *
     */
    ///@{

    MinHeap* nodes;                     // Mesh Node list
    Element **elements;               // Mesh Elements list
    Condition **dirichlet_conditions; // Mesh Dirichelet Conditions list
    Condition **neumann_conditions;   // Mesh Nueman Conditions list
    ///@}

public:
    Mesh() {}

    ~Mesh()
    {
        free(nodes);
        free(elements);
        free(dirichlet_conditions);
        free(neumann_conditions);
    }

    void set_problem_data(float k, float Q)
    {
        problem_data[THERMAL_CONDUCTIVITY] = k;
        problem_data[HEAT_SOURCE] = Q;
    }
    float get_problem_data(parameter position)
    {
        return problem_data[position];
    }

    void set_quantities(int num_nodes, int num_elements, int num_dirichlet, int num_neumann)
    {
        quantities[NUM_NODES] = num_nodes;
        quantities[NUM_ELEMENTS] = num_elements;
        quantities[NUM_DIRICHLET] = num_dirichlet;
        quantities[NUM_NEUMANN] = num_neumann;
    }
    int get_quantity(quantity position)
    {
        return quantities[position];
    }

    void init_arrays()
    {
        nodes = new MinHeap(quantities[NUM_NODES]);
        elements = (Element **)malloc(sizeof(Element *) * quantities[NUM_ELEMENTS]);
        dirichlet_conditions = (Condition **)malloc(sizeof(Condition *) * quantities[NUM_DIRICHLET]);
        neumann_conditions = (Condition **)malloc(sizeof(Condition *) * quantities[NUM_NEUMANN]);
    }

    void insert_node(Node *node, int position)
    {
        nodes->insert(node);
    }

    Node *get_node(int position)
    {
        Node * test = nodes->getNodeById(position);
        return test;
    }

    void zaxs(Element *element, int position)
    {
        elements[position] = element;
    }

    void insert_element(Element *element, int position)
    {
        elements[position] = element;
    }
    Element *get_element(int position)
    {
        return elements[position];
    }

    void insert_dirichlet_condition(Condition *dirichlet_condition, int position)
    {
        dirichlet_conditions[position] = dirichlet_condition;
    }

    Condition *get_dirichlet_condition(int position)
    {
        return dirichlet_conditions[position];
    }

    bool does_node_have_dirichlet_condition(int id)
    {
        bool ans = false;
        for (int i = 0; i < quantities[NUM_DIRICHLET]; i++)
            if (dirichlet_conditions[i]->get_node()->get_ID() == id)
            {
                ans = true;
                break;
            }
        return ans;
    }

    void insert_neumann_condition(Condition *neumann_condition, int position)
    {
        neumann_conditions[position] = neumann_condition;
    }

    Condition *get_neumann_condition(int position)
    {
        return neumann_conditions[position];
    }

    void report()
    {
        cout << "Problem Data\n**********************\n";
        cout << "Thermal Conductivity: " << problem_data[THERMAL_CONDUCTIVITY] << "\n";
        cout << "Heat Source: " << problem_data[HEAT_SOURCE] << "\n\n";
        cout << "Quantities\n***********************\n";
        cout << "Number of nodes: " << quantities[NUM_NODES] << "\n";
        cout << "Number of elements: " << quantities[NUM_ELEMENTS] << "\n";
        cout << "Number of dirichlet boundary conditions: " << quantities[NUM_DIRICHLET] << "\n";
        cout << "Number of neumann boundary conditions: " << quantities[NUM_NEUMANN] << "\n\n";
        cout << "List of nodes\n**********************\n";
        for (int i = 0; i < quantities[NUM_NODES]; i++)
            cout << "Node: " << nodes->getNodeById(i)->get_ID() << ", x= " <<  nodes->getNodeById(i)->get_x_coordinate() << ", y= " << nodes->getNodeById(i)->get_y_coordinate() << ", z= " <<  nodes->getNodeById(i)->get_z_coordinate() << "\n";
        cout << "\nList of elements\n**********************\n";
        for (int i = 0; i < quantities[NUM_ELEMENTS]; i++)
        {
            cout << "Element: " << elements[i]->get_ID() << ", Node 1= " << elements[i]->get_node1()->get_ID();
            cout << ", Node 2= " << elements[i]->get_node2()->get_ID() << ", Node 3= " << elements[i]->get_node3()->get_ID() << "Node 4= " << elements[i]->get_node4()->get_ID() << "\n";
        }
        cout << "\nList of Dirichlet boundary conditions\n**********************\n";
        for (int i = 0; i < quantities[NUM_DIRICHLET]; i++)
            cout << "Condition " << i + 1 << ": " << dirichlet_conditions[i]->get_node()->get_ID() << ", Value= " << dirichlet_conditions[i]->get_value() << "\n";
        cout << "\nList of Neumann boundary conditions\n**********************\n";
        for (int i = 0; i < quantities[NUM_NEUMANN]; i++)
            cout << "Condition " << i + 1 << ": " << neumann_conditions[i]->get_node()->get_ID() << ", Value= " << neumann_conditions[i]->get_value() << "\n";
        cout << "\n";
    }
};