#ifndef SSTREE_H
#define SSTREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include <limits>
#include <fstream>

#include "params.h"
#include "Point.h"


struct Pair {
    Point point;
    NType distance;

    Pair(const Point& p, NType d) : point(p), distance(d) {}
};

struct Comparator {
    bool operator()(const Pair& a, const Pair& b) const {
        return a.distance < b.distance; // max-heap basado en distancia
    }
};

class SsNode {
private:
    NType varianceAlongDirection(const std::vector<Point>& centroids, size_t direction) const;
    size_t minVarianceSplit(size_t coordinateIndex);
    
public:
    virtual ~SsNode() = default;

    Point centroid; 
    NType radius;
    SsNode* parent = nullptr;

    virtual bool isLeaf() const = 0;
    virtual std::vector<Point> getEntriesCentroids() const = 0;
    virtual void sortEntriesByCoordinate(size_t coordinateIndex) = 0;
    virtual std::pair<SsNode*, SsNode*> split() = 0;
    virtual bool intersectsPoint(const Point& point) const {
        return distance(this->centroid, point) <= this->radius;
    }

    virtual void updateBoundingEnvelope() = 0;
    size_t directionOfMaxVariance() const;
    size_t findSplitIndex();

    virtual pair<SsNode*,SsNode*> insert(const Point& point) = 0;

    bool test(bool isRoot = false) const;
    void print(size_t indent) const;

    virtual void FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const = 0;


    virtual void saveToStream(std::ostream &out) const = 0;
    virtual void loadFromStream(std::istream &in) = 0;
};

class SsInnerNode : public SsNode {
private:
    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    SsInnerNode() = default;
    SsInnerNode(size_t d);
    std::pair<SsNode*, SsNode*> split() override;
    std::vector<SsNode*> children;

    SsNode* findClosestChild(const Point& target) const;
    bool isLeaf() const override { return false; }
    void updateBoundingEnvelope() override;

    pair<SsNode*,SsNode*> insert(const Point& point) override;

    void FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const override;

    virtual void saveToStream(std::ostream &out) const override;
    virtual void loadFromStream(std::istream &in) override;
};

class SsLeaf : public SsNode {
private:

    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    SsLeaf() = default;
    SsLeaf(size_t d);
    std::pair<SsNode*, SsNode*> split() override;
    std::vector<std::string> paths;
    std::vector<Point> points;

    bool isLeaf() const override { return true; }
    void updateBoundingEnvelope() override;

    pair<SsNode*,SsNode*> insert(const Point& point) override;

    void FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const override;

    virtual void saveToStream(std::ostream &out) const override;
    virtual void loadFromStream(std::istream &in) override;
};


class SsTree {
private:
    SsNode* root;
    SsNode* search(SsNode* node, const Point& target);
    SsNode* searchParentLeaf(SsNode* node, const Point& target);

public:
    SsTree() : root(nullptr) {}
    ~SsTree() {
        delete root;
    }

    size_t D;

    NType nivel0() const {
        return root->radius;
    }

    NType nivel1() const {
        NType sum = 0;
        for (auto child : ((SsInnerNode*)root)->children) {
            sum += child->radius;
            cout << child->radius << endl;  
        }
        return sum;
    }

    void setD(size_t d) {
        D = d;
    }
    
    void insert(const Point& point);
    void insert( Point& point, const std::string& path);
    void build (const std::vector<Point>& points);
    std::vector<string> kNNQuery(const Point& center, size_t k) const;

    void print() const;
    void test() const;

    void saveToFile(const std::string &filename) const;
    void loadFromFile(const std::string &filename);
};

#endif // !SSTREE_H