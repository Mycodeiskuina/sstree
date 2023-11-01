#include "SStree.h"
const long long inf = 1e18;


NType SsNode::varianceAlongDirection(const std::vector<Point>& values, size_t direction) const {
    NType mean = 0;
    for (const Point& value : values) {
        mean += value[direction];
    }
    mean /= values.size();

    NType variance = 0;
    for (const Point& value : values) {
        variance += pow(value[direction] - mean, 2);
    }
    variance /= values.size();

    return variance;
}

size_t SsNode::minVarianceSplit(size_t coordinateIndex){
    size_t splitIndex = Settings::m;
    NType minVariance = inf;
    vector<Point> points = getEntriesCentroids();

    for (size_t i = Settings::m; i <= Settings::M - Settings::m; ++i) {
        std::vector<Point> left(points.begin(), points.begin() + i);
        std::vector<Point> right(points.begin() + i, points.end());

        NType leftVariance = varianceAlongDirection(left, coordinateIndex);
        NType rightVariance = varianceAlongDirection(right, coordinateIndex);
        NType variance = leftVariance + rightVariance;

        if (variance < minVariance) {
            minVariance = variance;
            splitIndex = i;
        }
    }

    return splitIndex;
}

size_t SsNode::directionOfMaxVariance() const { 
    std::vector<Point> centroids = getEntriesCentroids();

    size_t direction = 0;
    NType maxVariance = 0;
    for (size_t i = 0; i < centroids[0].dim(); ++i) {
        NType variance = varianceAlongDirection(centroids, i);
        if (variance > maxVariance) {
            maxVariance = variance;
            direction = i;
        }
    }

    return direction;
}

size_t SsNode::findSplitIndex() {
    size_t coordinateIndex = directionOfMaxVariance();
    sortEntriesByCoordinate(coordinateIndex);
    return minVarianceSplit(coordinateIndex);
}


std::vector<Point> SsInnerNode::getEntriesCentroids() const {
    std::vector<Point> centroids;
    for (const SsNode* child : children) {
        centroids.push_back(child->centroid);
    }
    return centroids;
}

void SsInnerNode::sortEntriesByCoordinate(size_t coordinateIndex) {
    std::sort(children.begin(), children.end(), [coordinateIndex](const SsNode* a, const SsNode* b) {
        return a->centroid[coordinateIndex] < b->centroid[coordinateIndex];
    });
}

SsNode* SsInnerNode::findClosestChild(const Point& target) const { 
    NType minDistance = inf;
    SsNode* closestChild = nullptr;
    for (SsNode* child : children) {
        NType distance = ::distance(child->centroid, target);
        if (distance < minDistance) {
            minDistance = distance;
            closestChild = child;
        }
    }
    return closestChild;
}

void SsInnerNode::updateBoundingEnvelope() {
    vector<Point> points = getEntriesCentroids();
    this->centroid = Point(points[0].dim());
    for (size_t i = 0; i < centroid.dim(); ++i) {
        for (const Point& point : points) {
            this->centroid[i] += point[i];
        }
    }
    
    this->centroid /= points.size();

    radius = 0;
    for (auto child : children) {
        NType distance = ::distance(this->centroid, child->centroid) + child->radius;
        if (distance > radius) {
            radius = distance;
        }
    }
}


std::vector<Point> SsLeaf::getEntriesCentroids() const {
    return points;
}

void SsLeaf::sortEntriesByCoordinate(size_t coordinateIndex) {
    std::sort(points.begin(), points.end(), [coordinateIndex](const Point& a, const Point& b) {
        return a[coordinateIndex] < b[coordinateIndex];
    });
}


void SsLeaf::updateBoundingEnvelope() { 
    vector<Point> points = getEntriesCentroids();
    this->centroid = Point(points[0].dim());
    for (size_t i = 0; i < centroid.dim(); ++i) {
        for (const Point& point : points) {
            this->centroid[i] += point[i];
        }
    }
    this->centroid /= points.size();

    radius = 0;
    for (const Point& point : points) {
        NType distance = ::distance(this->centroid, point);
        if (distance > radius) {
            radius = distance;
        }
    }
}

std::pair<SsNode*, SsNode*> SsLeaf::split() {
    size_t splitIndex = findSplitIndex();
    std::vector<Point> left(points.begin(), points.begin() + splitIndex + 1);
    std::vector<Point> right(points.begin() + splitIndex + 1, points.end());
    SsLeaf* leftNode = new SsLeaf();
    SsLeaf* rightNode = new SsLeaf();
    leftNode->points = left;
    rightNode->points = right;

    leftNode->updateBoundingEnvelope();
    rightNode->updateBoundingEnvelope();
    
    return std::make_pair(leftNode, rightNode);
}

std::pair<SsNode*, SsNode*> SsInnerNode::split() {
    size_t splitIndex = findSplitIndex();

    std::vector<SsNode*> left(children.begin(), children.begin() + splitIndex + 1);
    std::vector<SsNode*> right(children.begin() + splitIndex + 1, children.end());

    SsInnerNode* leftNode = new SsInnerNode();
    SsInnerNode* rightNode = new SsInnerNode();
    leftNode->children = left;
    rightNode->children = right;

    for (SsNode* child : left) {
        child->parent = leftNode;
    }

    for (SsNode* child : right) {
        child->parent = rightNode;
    }

    // Actualiza radio y centroide
    leftNode->updateBoundingEnvelope();
    rightNode->updateBoundingEnvelope();

    return std::make_pair(leftNode, rightNode);
}

pair<SsNode*,SsNode*> SsInnerNode::insert(const Point& point) {
    SsNode* closestChild = findClosestChild(point);
    pair<SsNode*,SsNode*> newChilds = closestChild->insert(point);
    pair<SsNode*, SsNode*> splitNodes;
    if (newChilds.first != nullptr) {
        children.erase(std::remove(children.begin(), children.end(), closestChild), children.end());
        children.push_back(newChilds.first);
        children.push_back(newChilds.second);
        if (children.size() > Settings::M) {
            splitNodes = split();
            splitNodes.first->parent = this;
            splitNodes.second->parent = this;
        }
    }
    closestChild->updateBoundingEnvelope();
    updateBoundingEnvelope();
    return splitNodes;
}

pair<SsNode*, SsNode*> SsLeaf::insert(const Point& point) {
    points.push_back(point);
    updateBoundingEnvelope();
    std::pair<SsNode*, SsNode*> splitNodes;
    if (points.size() > Settings::M) {
        splitNodes = split();
        splitNodes.first->parent = parent;
        splitNodes.second->parent = parent;
    }
    return splitNodes;
}


void SsTree::insert(const Point& point){
    if (!root) {
        root = new SsLeaf();
        dynamic_cast<SsLeaf*>(root)->points.push_back(point);
        root->updateBoundingEnvelope();
        root->parent = nullptr;
    }else{
        pair<SsNode*,SsNode*> newChilds = root->insert(point);
        if (newChilds.first != nullptr) {
            root = new SsInnerNode();
            dynamic_cast<SsInnerNode*>(root)->children.push_back(newChilds.first);
            newChilds.first->parent = root;
            dynamic_cast<SsInnerNode*>(root)->children.push_back(newChilds.second);
            newChilds.second->parent = root;
            dynamic_cast<SsInnerNode*>(root)->updateBoundingEnvelope();
            dynamic_cast<SsInnerNode*>(root)->parent = nullptr;
        }
    }
    //print();
    //cout<<"----------------------------------------------------------------------------"<<endl;
}

void SsTree::insert(Point& point, const std::string& path){
    point.path = "../" + path;
    if (!root) {
        root = new SsLeaf();
        dynamic_cast<SsLeaf*>(root)->points.push_back(point);
        root->updateBoundingEnvelope();
        root->parent = nullptr;
    }else{
        pair<SsNode*,SsNode*> newChilds = root->insert(point);
        if (newChilds.first != nullptr) {
            root = new SsInnerNode();
            dynamic_cast<SsInnerNode*>(root)->children.push_back(newChilds.first);
            newChilds.first->parent = root;
            dynamic_cast<SsInnerNode*>(root)->children.push_back(newChilds.second);
            newChilds.second->parent = root;
            dynamic_cast<SsInnerNode*>(root)->updateBoundingEnvelope();
            dynamic_cast<SsInnerNode*>(root)->parent = nullptr;
        }
    }
}

SsNode* SsTree::search(SsNode* node, const Point& target){
    if (node->isLeaf()) {
        return node;
    } else {
        SsInnerNode* inner = dynamic_cast<SsInnerNode*>(node);
        SsNode* closestChild = inner->findClosestChild(target);
        return search(closestChild, target);
    }
}

SsNode* SsTree::searchParentLeaf(SsNode* node, const Point& target) {
    if (node->isLeaf()) {
        return node;
    } else {
        SsInnerNode* inner = dynamic_cast<SsInnerNode*>(node);
        SsNode* closestChild = inner->findClosestChild(target);
        return searchParentLeaf(closestChild, target);
    }
}




void SsTree::build(const std::vector<Point>& points){
    for (const Point& point : points) {
        insert(point);
    }
}


void SsLeaf::FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const{
    Point M = this->centroid;
    for (const Point& point : points) {
        if (distance(M,q) - distance(M,point) > Dk) {
            continue;
        } else if (distance(M,point) - distance(M,q) > Dk) {
            continue;
        } else {
            if (distance(point, q) < Dk) {
                if (L.size() == k) {
                    L.pop();
                }
                Pair p(point, distance(point, q));
                L.push(p);
                if (L.size() == k) {
                    Dk = L.top().distance;
                }
            }
        }
    }
}


void SsInnerNode::FNDFTrav(const Point& q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator>& L, NType& Dk) const{
    for (const SsNode* child : children) { 
        child->FNDFTrav(q, k, L, Dk);
    }
}


vector<string> SsTree::kNNQuery(const Point& center, size_t k) const{
    std::priority_queue<Pair, std::vector<Pair>, Comparator> L;
    NType Dk = inf;
    root->FNDFTrav(center, k, L, Dk);
    vector<string> paths;
    while (!L.empty()) {
        paths.push_back(L.top().point.path);
        L.pop();
    }
    return paths;
}

bool SsNode::test(bool isRoot) const {
    size_t count = 0;
    if (this->isLeaf()) {
        const SsLeaf* leaf = dynamic_cast<const SsLeaf*>(this);
        count = leaf->points.size();

        for (const Point& point : leaf->points) {
            if (distance(this->centroid, point) > this->radius) {
                std::cout << "Point outside node radius detected." << std::endl;
                return false;
            }
        }
    } else {
        const SsInnerNode* inner = dynamic_cast<const SsInnerNode*>(this);
        count = inner->children.size();

        for (const SsNode* child : inner->children) {
            if (distance(this->centroid, child->centroid) > this->radius) {
                std::cout << "Child centroid outside parent radius detected." << std::endl;
                return false;
            }
            if (!child->test(false)) {
                return false;
            }
        }
    }

    if (!isRoot && (count < Settings::m || count > Settings::M)) {
        std::cout << "Invalid number of children/points detected." << std::endl;
        return false;
    }

    if (!isRoot && !parent) {
        std::cout << "Node without parent detected." << std::endl;
        return false;
    }

    return true;
}

void SsTree::test() const {
    bool result = root->test(true);

    if (root->parent) {
        std::cout << "Root node parent pointer is not null!" << std::endl;
        result = false;
    }

    if (result) {
        std::cout << "SS-Tree is valid!" << std::endl;
    } else {
        std::cout << "SS-Tree has issues!" << std::endl;
    }
}


void SsNode::print(size_t indent = 0) const {
    for (size_t i = 0; i < indent; ++i) {
        std::cout << "  ";
    }

    // Imprime información del nodo.
    std::cout << "Centroid: " << centroid << ", Radius: " << radius;
    if (isLeaf()) {
        const SsLeaf* leaf = dynamic_cast<const SsLeaf*>(this);
        std::cout << ", Points: [ ";
        for (const Point& p : leaf->points) {
            std::cout << p << " ";
        }
        std::cout << "]";
    } else {
        std::cout << std::endl;
        const SsInnerNode* inner = dynamic_cast<const SsInnerNode*>(this);
        for (const SsNode* child : inner->children) {
            child->print(indent + 1); 
        }
    }
    std::cout << std::endl;
}
void SsTree::print() const {
    if (root) {
        root->print();
    } else {
        std::cout << "Empty tree." << std::endl;
    }
}


void SsLeaf::saveToStream(std::ostream &out) const {
    //cout << "saveToStream Leaf" << endl;
    // Guardar centroid
    centroid.saveToFile(out, Settings::D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char*>(&radius_), sizeof(radius_));

    // Guardar el numero de puntos
    size_t numPoints = points.size();
    out.write(reinterpret_cast<const char*>(&numPoints), sizeof(numPoints));

    // Guardar los puntos
    for (const auto& point : points) {
        point.saveToFile(out,  Settings::D);
    }

    // Guardar las rutas (paths)
    size_t numPaths = points.size();
    out.write(reinterpret_cast<const char*>(&numPaths), sizeof(numPaths));
    for (const auto& p : points) {
        size_t pathLength = p.path.size();
        out.write(reinterpret_cast<const char*>(&pathLength), sizeof(pathLength));
        out.write(p.path.c_str(), (long) pathLength);
        //cout<<"path: "<<p.path<<endl;
    }
}

void SsInnerNode::saveToStream(std::ostream &out) const {
    //cout<<"saveToStream Inner"<<endl;
    // Guardar centroid
    centroid.saveToFile(out,  Settings::D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char*>(&radius_), sizeof(radius_));

    // Guardar si apunta a nodos hoja
    bool pointsToLeafs = children[0]->isLeaf();
    out.write(reinterpret_cast<const char*>(&pointsToLeafs), sizeof(pointsToLeafs));

    // Guardar la cantidad de hijos para saber cuántos nodos leer después
    size_t numChildren = children.size();
    out.write(reinterpret_cast<const char*>(&numChildren), sizeof(numChildren));

    // Guardar los hijos
    for (const auto& child : children) {
        child->saveToStream(out);
    }
}

void SsInnerNode::loadFromStream(std::istream &in) {
    // Leer centroid
    centroid.readFromFile(in,  Settings::D);

    // leer el valor del radio
    float radius_ = 0;
    in.read(reinterpret_cast<char*>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // leer si apunta a hojas o nodos internos
    bool pointsToLeaf = false;
    in.read(reinterpret_cast<char*>(&pointsToLeaf), sizeof(pointsToLeaf));

    // leer cantidad de hijos
    size_t numChildren;
    in.read(reinterpret_cast<char*>(&numChildren), sizeof(numChildren));

    // leer hijos
    for (size_t i = 0; i < numChildren; ++i) {
        SsNode* child = pointsToLeaf ? static_cast<SsNode*>(new SsLeaf()) : static_cast<SsNode*>(new SsInnerNode());
        child->loadFromStream(in);
        children.push_back(child);
    }
}

void SsLeaf::loadFromStream(std::istream &in) {
    //cout<<"loadFromStream Leaf"<<endl;
    // Leer centroid
    centroid.readFromFile(in,  Settings::D);

    // Leer radio
    float radius_ = 0;
    in.read(reinterpret_cast<char*>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // Leer numero de puntos
    size_t numPoints;
    in.read(reinterpret_cast<char*>(&numPoints), sizeof(numPoints));

    // Leer puntos
    points.resize(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        points[i].readFromFile(in,  Settings::D);
    }

    // Leer rutas (paths)
    size_t numPaths;
    in.read(reinterpret_cast<char*>(&numPaths), sizeof(numPaths));
    for (size_t i = 0; i < numPaths; ++i) {
        size_t pathLength;
        in.read(reinterpret_cast<char*>(&pathLength), sizeof(pathLength));
        char* buffer = new char[pathLength + 1];
        in.read(buffer, (long) pathLength);
        buffer[pathLength] = '\0';
        points[i].path = std::string(buffer);
        //cout << "path: " << points[i].path << endl;
        delete[] buffer;
    }
}

void SsTree::saveToFile(const std::string &filename) const {
    //cout<<"saveToFile"<<endl;
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Cannot open file for writing");
    }

    // Guardar las dimensiones de la estructura
    out.write(reinterpret_cast<const char*>(&D), sizeof(D));

    // Guardar si el root es hija o nodo interno
    bool isLeaf = root->isLeaf();
    out.write(reinterpret_cast<const char*>(&isLeaf), sizeof(isLeaf));

    // Guardar el resto de la estructura
    root->saveToStream(out);
    out.close();
}

void SsTree::loadFromFile(const std::string &filename) {
    //cout<<"loadFromFile"<<endl;
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open file for reading");
    }
    if (root) {
        delete root;
        root = nullptr;
    }

    // Aquí se asume que el primer valor determina las dimensiones
    in.read(reinterpret_cast<char*>(&D), sizeof(D));

    // El segundo valor determina si el root es hoja
    bool isLeaf;
    in.read(reinterpret_cast<char*>(&isLeaf), sizeof(isLeaf));
    if (isLeaf) {
        root = new SsLeaf();
    } else {
        root = new SsInnerNode();
    }
    root->loadFromStream(in);
    in.close();
}

