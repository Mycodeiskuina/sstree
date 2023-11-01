#include <iostream>
#include <vector>
#include <random>
#include "Point.h"
#include "SStree.h"
#include <hdf5/serial/H5Cpp.h>
#include <nlohmann/json.hpp>
#include <fstream>

string myR = "/mnt/c/labo-eda/code/";

struct ImageData {
    Point embedding;
    std::string path;
};


std::vector<ImageData> readEmbeddingsFromJson(const std::string& FILE_NAME) {
    std::vector<ImageData> data;

    try {
        std::ifstream file(FILE_NAME);
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open JSON file.");
        }

        nlohmann::json jsonData;
        file >> jsonData;
        std::vector<std::vector<float>> features = jsonData["features"].template get<std::vector<std::vector<float>>>();
        std::vector<std::string> paths = jsonData["paths"].template get<std::vector<std::string>>();

        for (size_t i = 0; i < paths.size(); ++i) {
            Point embedding(features[i].size());
            for (size_t j = 0; j < features[i].size(); ++j) {
                embedding[j] = features[i][j];
            }
            data.push_back({embedding, paths[i]});
        }
        
        file.close();
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return data;
}

std::vector<ImageData> readEmbeddingsFromHDF5(const H5std_string& FILE_NAME) {
    std::vector<ImageData> data;
    try{
        H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);

        // Leer dataset "features"
        H5::DataSet featuresDataset = file.openDataSet("features");
        H5::DataSpace featuresSpace = featuresDataset.getSpace();
        const int num_embeddings = featuresSpace.getSimpleExtentNpoints();
        double* embeddingData = new double[num_embeddings];
        featuresDataset.read(embeddingData, H5::PredType::NATIVE_DOUBLE);
        // Leer dataset "paths"
        H5::DataSet pathsDataset = file.openDataSet("paths");
        H5::StrType strType = pathsDataset.getStrType();
        H5::DataSpace pathsSpace = pathsDataset.getSpace();
        const int num_paths = pathsSpace.getSimpleExtentNpoints();
        std::vector<H5std_string> paths(num_paths);
        pathsDataset.read(paths.data(), strType);
        for (int i = 0; i < num_paths; ++i) {
            Point embedding(num_embeddings);
            for (int j = 0; j < num_embeddings; ++j) {
                embedding[j] = embeddingData[i * num_embeddings + j];
            }
            data.push_back({embedding, paths[i]});
        }
        delete[] embeddingData;
        file.close();
    } catch(H5::Exception& error) {
        std::cerr << error.getCDetailMsg() << std::endl;
    }

    return data;
}

int main() {
    const std::string FILE_NAME = "../embedding.json";
    std::vector<ImageData> data = readEmbeddingsFromJson(FILE_NAME);
    SsTree tree;
    for (ImageData& item : data) {
        item.embedding;
        tree.insert(item.embedding, item.path);
    }
    cout<<tree.nivel0()<<endl;
    cout<<tree.nivel1()<<endl;
    tree.test();
    std::string filename = "../embbeding.dat";
    tree.saveToFile(filename);
}
