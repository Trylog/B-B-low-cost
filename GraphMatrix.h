//
// Created by micha on 08.05.2023.
//

#include "Edge.h"

namespace std {

    class GraphMatrix {
    public:
        int v{}; //liczba wierzchołków
        int** matrix;
        explicit GraphMatrix(int v);
        explicit GraphMatrix(string path);
        ~GraphMatrix();
        void insert(Edge edge);
        void display();
    };

} // std

