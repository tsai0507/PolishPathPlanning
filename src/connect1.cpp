#include <iostream>
#include <X11/Xlib.h>
#include <pcl/console/print.h>
#include <pcl/console/parse.h>
#include "Path_Generate_Algorithm.h"

int main(int argc, char **argv)
{
    //deal with thread problem of pcl occured sometime
    XInitThreads();
    if (argc < 2)
    {
        cout << "./slicing_method cad_name.pcd" << endl;
        return (-1);
    }
    std::vector<int> pcd_file_indices = pcl::console::parse_file_extension_argument(argc, argv, ".pcd");
    if (pcd_file_indices.size() != 1)
    {
        cout << "./slicing_method cad_name.pcd" << endl;
        return (-1);
    }

    std::string configFile = "../config.txt";
    path_generater path_planner = {configFile, argv[pcd_file_indices[0]]};
    path_planner.GenPath();
    path_planner.getPath();
    path_planner.show();

    return 0;
}