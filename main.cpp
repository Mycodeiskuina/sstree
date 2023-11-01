#include <iostream>
#include <vector>
#include <random>
#include "SStree.h"

int main() {
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-10.0, 10.0);

    std::vector<Point> points(500, Point(50));
    for (int i = 0; i < 500; i++) {
        for (int j = 0; j < 50; j++) {
            points[i][j] = dis(gen);
        }
    }

    // Insert the points into the SStree
    SsTree tree;
    for (auto& point : points) {
        tree.insert(point);
    }
    //tree.print();
    cout << tree.nivel0() << endl;
    cout << tree.nivel1() << endl;
    tree.test();
    //std::string filename = "sstree.dat";
    //tree.saveToFile(filename);

    tree = SsTree();    // Clean the tree
    //tree.loadFromFile(filename);
    
    return 0;
}
