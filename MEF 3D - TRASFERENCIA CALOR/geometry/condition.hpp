/**
 * @file geometry/condition.hpp
 *
 * @author Enmanuel Amaya, MSc. 2023 *
 * @author Modified By: Grupo 09 - sec 01
 *
 * @brief Condition Type
 * @version 1
 * @date 2023-06-30
 *
 * A condition is a relation of a node and boundary condition value 
 * has a setter and getter for the private objects
 * 
 * Specifies which nodes should be considered to apply the Dirichlet and Neumman boundary conditions. */
class Condition {
    private:
        Node* node;
        float value;

    public:
        Condition(Node* node_to_assign, float value_to_assign){
            node = node_to_assign;
            value = value_to_assign;
        }

        void set_node(Node* node_to_assign){
            node = node_to_assign;
        }
        
        Node* get_node(){
            return node;
        }

        void set_value(float value_to_assign){
            value = value_to_assign;
        }
        float get_value(){
            return value;
        }
};
