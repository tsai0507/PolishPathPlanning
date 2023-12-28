#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <Eigen/Core>

void write_txt(std::string path_name, std::vector<Eigen::VectorXd> data);

int main(int argc, char **argv)
{
   //open .txt file which save path
    std::string line, TxtFileName;
    std::vector<Eigen::VectorXd> KeyCartesianPoint, CartesianPoint;
    std::ifstream inputFile("../ft_test.txt");

    if (!inputFile.is_open()) {
        std::cout << "Can't open file:" << TxtFileName << std::endl;
        return false;
    }
    while (std::getline(inputFile, line)) 
    {
        double value;
        std::istringstream iss(line);
        Eigen::VectorXd PointC_Eigen(6);

        //read Cartesian data(x,y,z,rx,ry,rz) from .txt
        for(size_t index = 0; index < 6; ++index){
            iss >> value;
            PointC_Eigen[index] = value;
        } 
        //convert Cartesian space data to joint space data
        KeyCartesianPoint.push_back(PointC_Eigen);
    } inputFile.close();
    

    

    return 0;
}

void write_txt(std::string path_name, std::vector<Eigen::VectorXd> data)
{
    size_t datalen = data.size();
    std::ofstream outputFile;
    outputFile.open(path_name);
    if (outputFile.is_open()){
        for (size_t index = 0; index < datalen; index++){
            for (size_t j = 0; j < 6; j++)
                outputFile << data[index][j] << " ";
            outputFile << std::endl;
        }
        outputFile.close();
        std::cout << "Output File saved: " << path_name << std::endl;
    }
    else
        std::cerr << "Unable to open file: " << path_name << std::endl;

}