/**
 * @name 3D MEF CHANGE #2
 *
 * - In 3D MESH each elements must have 4 nodes 
 */
/**
 * @brief MESH Element
 *
 * A colection of nodes witch represent a single element in MESH
 * - 1D MESH : 2 NODES
 * - 2D MESH : 3 NODES
 * - 3D MESH : 4 NODES
 *
 */
class Element
{
private:
   int ID;
    Node *node1;
    Node *node2;
    Node *node3;
    Node *node4;

public:
    /**
     * @brief Element constructor
     *
     * @param identifier
     * @param first_node
     * @param second_node
     * @param third_node
     * @param fourth_node
     */
    Element(int identifier, Node *first_node, Node *second_node, Node *third_node, Node *fourth_node)
    {
        ID = identifier;
        node1 = first_node;
        node2 = second_node;
        node3 = third_node;
        node4 = fourth_node;
    }

    /**
     * @name Basic setters and getters
     */
    
    void set_ID(int identifier)
    {
        ID = identifier;
    }
   int get_ID()
    {
        return ID;
    }

    void set_node1(Node *node)
    {
        node1 = node;
    }
    Node *get_node1()
    {
        return node1;
    }

    void set_node2(Node *node)
    {
        node2 = node;
    }
    Node *get_node2()
    {
        return node2;
    }

    void set_node3(Node *node)
    {
        node3 = node;
    }
    Node *get_node3()
    {
        return node3;
    }

    void set_node4(Node *node)
    {
        node4 = node;
    }
    Node *get_node4()
    {
        return node4;
    }
};