#ifndef PATH_GENERATION
#define PATH_GENERATION

#include <time.h>
#include <string>
#include <vector>
#include <map>

// #include <Eigen/Dense>
// #include <Eigen/Eigenvalues>

#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/PointIndices.h>
#include <pcl/common/centroid.h>
#include <pcl/common/transforms.h>
#include <pcl/common/common.h>
#include <pcl/common/impl/angles.hpp>
#include <pcl/features/normal_3d.h>
#include <pcl/search/kdtree.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/principal_curvatures.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/surface/mls.h>

#include "Spline.h"

typedef pcl::PointCloud<pcl::PointXYZRGB> PointCloudType;
typedef pcl::PointXYZRGB PointType;

class path_generater{
public:
    path_generater(){};
    path_generater(std::string cloud_name, double Radius);
    ~path_generater(){};

    void show();
    void voxel_down(const float x, const float y, const float z);
    void trans2center();
    void slicing_method();
    std::map<double, std::vector<double>>  
        insert_point(std::vector<int> indices, Eigen::Vector3f PlanePoint);
    void Set_kdtree();
    void estimate_normal();
    void smooth();
    void Contact_Path_Generation();
    std::vector<int> rangedX_index(int position);
    void get_coverage();

private:
    void path_track(Eigen::Vector3f plane_point);
    void dynamic_adjust_path(Spline* origin_path, Spline pre_path);
    int  compute_boundary(Spline path, Spline* boundary);
    void drawpath(Spline path,int r, int g, int b);
    //flag(1 draw)  key(0/1 up/down) 
    Eigen::Vector3f Area2Cloud(Eigen::Vector3d PathNode, bool flag, bool key);
    Eigen::Matrix4f compute_transform(pcl::PointXYZRGB point, float *principle_curvature);
    Eigen::Vector3d bisection(Eigen::Vector3d point, Spline* boundary, int itr);
    void compute_coverage(Eigen::Vector3d node, double radius);


    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, other_cloud;
    pcl::KdTreeFLANN<pcl::PointXYZRGB> kdtree;
    pcl::PointCloud<pcl::Normal>::Ptr cloud_with_normals;
    Eigen::Matrix4f transform; //camera to origin
    Eigen::Vector3f normal;
    double toolRadius, depth = 0.005, Adjust_Threshold = 1, toolthickness = 10;
    std::vector<Spline> Path_set;
    std::string file_name;
    std::vector<int> coverage_flag;
    int cloud_size;
};      


#endif