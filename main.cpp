///Michal Bernacki-Janson
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <utility>
#include <cstring>
#include <cstdint>
#include <chrono>
#include "GraphMatrix.h"

using namespace std;
int N = 0;
int *minSeq;

struct Node{
    vector<pair<int, int>> path;
    int** reducedMatrix;
    int bound{}, vertex{}, depth{};
};

Node* newNode(int** parentMatrix, vector<pair<int, int>> const& path, int depth, int i, int j){
    Node* node = new Node;
    node->path = path;
    if (depth != 0) node->path.emplace_back(i, j);

    node->reducedMatrix = new int*[N];
    for (int k = 0; k < N; ++k) {
        node->reducedMatrix[k] = new int[N];
        memcpy(node->reducedMatrix[k], parentMatrix[k],(N* sizeof(int)));
    }

    for (int k = 0; depth != 0 && k < N; k++){
        node->reducedMatrix[i][k] = INT32_MAX;
        node->reducedMatrix[k][j] = INT32_MAX;
    }

    node->reducedMatrix[j][0] = INT32_MAX;
    node->depth = depth;
    node->vertex = j;

    return node;
}

int calcCost(int** reducedMatrix){
    int cost = 0;
    int row[N];
    int col[N];

    fill_n(row, N, INT32_MAX);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (reducedMatrix[i][j] < row[i]) {
                row[i] = reducedMatrix[i][j];
            }
        }
    }
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (reducedMatrix[i][j] != INT32_MAX && row[i] != INT32_MAX) {
                reducedMatrix[i][j] -= row[i];
            }
        }
    }

    fill_n(col, N, INT32_MAX);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (reducedMatrix[i][j] < col[j]) col[j] = reducedMatrix[i][j];
        }
    }
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (reducedMatrix[i][j] != INT32_MAX && col[j] != INT32_MAX) reducedMatrix[i][j] -= col[j];
        }
    }

    for (int i = 0; i < N; i++){
        cost += (row[i] != INT32_MAX) ? row[i] : 0,
        cost += (col[i] != INT32_MAX) ? col[i] : 0;
    }
    return cost;
}

void printPath(vector<pair<int, int>> const& list){
    for (int i = 0; i < N; ++i) {
        minSeq[i] = list[i].first;
    }
}
///funckja porównywania w kolejce priorytetowej (lambda)
auto comp = [](const Node* n1, const Node* n2){return n1->bound > n2->bound;};

int solve(int** costMatrix){
    priority_queue<Node*, vector<Node*>, decltype(comp)> queue(comp);
    vector<pair<int, int>> v;
    Node* root = newNode(costMatrix, v, 0, -1, 0);

    root->bound = calcCost(root->reducedMatrix);
    queue.push(root);

    while (!queue.empty()){
        Node* min = queue.top();
        queue.pop();

        int i = min->vertex;
        if (min->depth == N - 1){
            min->path.emplace_back(i, 0);
            printPath(min->path);
            return min->bound;
        }

        for (int j = 0; j < N; j++){
            if (min->reducedMatrix[i][j] != INT32_MAX){
                Node* child = newNode(min->reducedMatrix, min->path,min->depth + 1, i, j);
                child->bound = min->bound + min->reducedMatrix[i][j] + calcCost(child->reducedMatrix);
                queue.push(child);
            }
        }
        delete[] min->reducedMatrix;
        delete min;
    }
}

int main()
{
    minSeq = new int[40];
    fstream in;
    fstream out;
    //
    ///otwieranie plik konfiguracyjnego, wczytawnie kolejnych lini, obliczanie ich liczby i wykrywanie nazwy pliku wyjściowego
    in.open("start.ini", ios::in);
    if(!in.is_open()){
        cout<<"Nie można otworzyć pliku konfiguracyjnego!"<<endl;
        return 1;
    }
    int tests = 0;
    string line;
    string csvName;
    vector<string> namesOfTests;
    while(getline(in, line)){
        if(in.eof()){
            csvName = line;
        } else{
            namesOfTests.push_back(line);
            tests++;
        }
    }
    in.close();
    //
    out.open(csvName, ios::out);
    if(!out.is_open()){
        cout<<"Nie można utworzyć pliku wynikowego"<<endl;
        return 0;
    }
    for(int t =0; t<tests;t++)
    {
        auto temp = namesOfTests.at(t).find(' ');
        //
        ///ekstrakcja nazwy pliku testu z lini
        GraphMatrix matrix = GraphMatrix(namesOfTests.at(t).substr(0, temp));
        //
        double timeAvg = 0.0;
        //
        ///ekstrakcja liczby powtórzeń danego testu
        int b1 = stoi(namesOfTests.at(t).substr(temp, namesOfTests.at(t).substr(temp).find(' ')-temp));
        //
        int minSum;

        for (int k = 0; k < b1; k++) {

            //----------------//
            ///główny algorytm

            auto t1 = std::chrono::high_resolution_clock::now();//pomiar czasu
            N = matrix.v;
            minSum = solve(matrix.matrix);

            //
            auto t2 = std::chrono::high_resolution_clock::now();
            auto tDiv = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
            //matrix.display();
            cout << endl << minSum << endl;
            for (int i = 0; i < matrix.v; ++i) {
                //cout << minSeq[i] << " ";
            }
            cout << tDiv.count();
            out << tDiv.count() << endl;
            timeAvg += tDiv.count();

            //if(k!=b1-1)delete[] minSeq;
            //----------------//
        }
        timeAvg /= (double) b1;
        out << namesOfTests.at(t).substr(0, temp) <<";" << timeAvg<<";" << minSum<< "; [";
        for (int i = 0; i < matrix.v-1; ++i) {
            out << minSeq[i] << " ";
        }
        out<<"0]"<<endl<<endl;

    }
    delete[] minSeq;
    out.close();
    //system("pause");
    return 0;
}
