
#include <iostream>
#include <fstream>
#include "GraphMatrix.h"

namespace std {
    GraphMatrix::GraphMatrix(int v) {
        this->v=v;
        matrix = new int*[v];
        for (int i = 0; i < v; ++i) {
            matrix[i] = new int[v];
            for (int j = 0; j < v; ++j) matrix[i][j]=0;
        }
    }

    GraphMatrix::GraphMatrix(string path) {
        fstream input;
        input.open(path, ios::in);
        if (!input.good()) cout<<"błąd pliku"<<endl;
        input >> v;
        matrix = new int*[v];

        for (int i = 0; i < v; ++i) {

            matrix[i] = new int[v];

            for (int j = 0; j < v; ++j) {
                int temp;
                input >> temp;
                matrix[i][j] = temp;
            }
        }
        input.close();
    }

    GraphMatrix::~GraphMatrix() {
        for (int i = 0; i < v; ++i) {
            delete[] matrix[i];
        }
        delete[] matrix;
    }

    void GraphMatrix::insert(Edge edge) {
        matrix[edge.n1][edge.n2]=edge.weight;
    }

    void GraphMatrix::display() {
        for (int i = 0; i < v; ++i) {
            for (int j = 0; j < v; ++j) {
                cout<<matrix[i][j]<<"; ";
            }
            cout<<endl;
        }
    }
} // std