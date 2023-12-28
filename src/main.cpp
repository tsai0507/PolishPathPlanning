#include <iostream>
#include "Path_Generate.h"
#include <pcl/console/print.h>
#include <pcl/console/parse.h>

int main(int argc, char **argv)
{
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

    double step_size = 0;
    std::cout << "Please enter the Radius: ";
    std::cin >> step_size;
    path_generater path_planner(argv[pcd_file_indices[0]], step_size);
    // path_planner.voxel_down(0.1, 1, 1);    
    // path_planner.trans2center();

    // path_planner.slicing_method();
    // path_planner.smooth();
    path_planner.estimate_normal();
    path_planner.Contact_Path_Generation();
    path_planner.get_coverage();
    path_planner.show();

    return 0;
}