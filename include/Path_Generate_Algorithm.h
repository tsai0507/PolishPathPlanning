#ifndef PATH_CONNECT
#define PATH_CONNECT

#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <chrono>
#include <memory>
#include <thread>
#include <functional>

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

// #define HANDEYEx -0.867557
// #define HANDEYEy 0.069088
// #define HANDEYEz 0.664201
// #define HANDEYErx -3.1270175
// #define HANDEYEry -0.040124
// #define HANDEYErz -1.6063578

#define HANDEYEx -0.764091
#define HANDEYEy 0.025886
#define HANDEYEz 0.663790
#define HANDEYErx -3.1270175
#define HANDEYEry -0.040124
#define HANDEYErz -1.6063578


typedef pcl::PointCloud<pcl::PointXYZRGB> PointCloudType;
typedef pcl::PointXYZRGB PointType;
typedef std::map<double, std::vector<double>> MAP;


// Base class which can genrate eual spacing path and output robotic path
class SectPath
{
public:
    SectPath(){};
    SectPath(std::string configName, std::string CloudFileName);

    void show();
    void estimate_normal();
    virtual void GenPath();
    void getPath();

protected:
    virtual void read_config(std::string filename);
    void trans2center();
    void remove_outlier();
    void smooth();
    void drawpath(Spline path, int r, int g, int b);

    //Path planning Alg.
    std::vector<int> rangedX_index(int position);
    MAP insert_point(std::vector<int> indices, Eigen::Vector3f PlanePoint);
    Spline OnePath(Eigen::Vector3f plane_point);
    
    //final path adjustment
    void HandEyeTransform(Eigen::VectorXf &point);
    void reduceRPY(std::vector<int> Index);
    void postion_smooth();
    void TransFlangeposition();

    // pcl using object
    PointCloudType::Ptr cloud, other_cloud;
    pcl::PointCloud<pcl::Normal>::Ptr normal_cloud;
    pcl::KdTreeFLANN<PointType> kdtree;

    Eigen::Matrix4f TransAlign; // camera to origin
    double toolRadius, PathResolution, RPYres;
    float EElen;
    bool ifAlign, ifSmooth, ifChangeRange, ifRemove;
    std::vector<Spline> Path_set;
    std::string cloud_name, pathFile;
    std::vector<Eigen::VectorXf> WayPointsList;
};


// derived class which can make dynamic adjustment
class path_generater : public SectPath 
{
public:
    path_generater(){};
    path_generater(std::string configName, std::string CloudFileName);
    void GenPath() override;

    typedef enum {left, right} Dir;
    struct thread_wrap_data{
        Dir key;
        int min_pt, max_pt;
        Spline CenPath;
        thread_wrap_data(){};
        thread_wrap_data(Dir _key, int _min_pt, int _max_pt, Spline _CenPath):
            key(_key), min_pt(_min_pt), max_pt(_max_pt) {CenPath = _CenPath;}
    }; typedef struct thread_wrap_data thread_wrap_data;

private:
    void read_config (std::string filename) override;

    // Alg dynamic adjustment
    Eigen::Matrix4f compute_transform(PointType point, float *principle_curvature);
    Eigen::Vector3f Area2Cloud(Eigen::Vector3d point, bool flag, Dir key);
    int compute_boundary(Spline path, std::shared_ptr<Spline>& boundary, Dir key); // use key to chose up/down
    Eigen::Vector3d bisection(Eigen::Vector3d PathNode, std::shared_ptr<Spline> boundary, int itr, Dir key);
    void dynamic_adjust_path(Spline *origin_path, std::shared_ptr<Spline> pre_path, Dir key);
    void inital_center_path(float position_X);
    // use for thread
    void thread_worker(thread_wrap_data data, std::vector<Spline> &_path_set);
   
    std::mutex mtx;
    double depth, Adjust_Threshold, toolthickness;
    bool Adjust;
};


#endif