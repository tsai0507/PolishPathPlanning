#ifndef PATH_CONNECT
#define PATH_CONNECT

#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>

// #include <Eigen/Dense>
// #include <Eigen/Eigenvalues>
// #include <Eigen/Geometry>
// pcl library contain Eigen

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

// #define HANDEYEx 0.792078
// #define HANDEYEy -0.042662
// #define HANDEYEz 0.6656017
// #define HANDEYErx -3.1531625
// #define HANDEYEry -0.048573
// #define HANDEYErz 1.609157

// #define HANDEYEx 0.17205
// #define HANDEYEy -0.847092
// #define HANDEYEz 0.664335
// #define HANDEYErx 3.1372645
// #define HANDEYEry -0.0359714
// #define HANDEYErz -3.0986405

// #define HANDEYEx -1.021639
// #define HANDEYEy 0.0538
// #define HANDEYEz 0.655478
// #define HANDEYErx 3.1105427
// #define HANDEYEry -0.0742821
// #define HANDEYErz -1.6327766

#define HANDEYEx -0.858533
#define HANDEYEy 0.075348
#define HANDEYEz 0.672533
#define HANDEYErx -3.138775
#define HANDEYEry -0.0405313
#define HANDEYErz -1.5707969

typedef pcl::PointCloud<pcl::PointXYZRGB> PointCloudType;
typedef pcl::PointXYZRGB PointType;
typedef std::map<double, std::vector<double>> MAP;

class SectPath
{
public:
    SectPath(){};
    SectPath(std::string configName, std::string CloudFileName, double Radius);

    void show();
    void estimate_normal();
    void GenPath();
    void getPath();

private:
    void read_config(std::string filename);
    void trans2center();
    void remove_outlier();
    void smooth();
    void drawpath(Spline path, int r, int g, int b);
    //Path planning Alg.
    std::vector<int> rangedX_index(int position);
    MAP insert_point(std::vector<int> indices, Eigen::Vector3f PlanePoint);
    void OnePath(Eigen::Vector3f plane_point);
    //final path adjustment
    void HandEyeTransform(Eigen::VectorXf &point);
    void reduceRPY(std::vector<int> Index);
    void TransFlangeposition();

    PointCloudType::Ptr cloud, other_cloud;
    pcl::PointCloud<pcl::Normal>::Ptr normal_cloud;
    pcl::KdTreeFLANN<PointType> kdtree;

    Eigen::Matrix4f TransAlign; // camera to origin
    double toolRadius, depth, Adjust_Threshold, toolthickness, PathResolution, RPYres;
    float EElen;
    bool ifAlign, ifSmooth, ifChangeRange, ifRemove;
    std::vector<Spline> Path_set;
    std::string cloud_name, pathFile;
    std::vector<int> coverage_flag;
    std::vector<Eigen::VectorXf> WayPointsList;
};


#endif