#include <iostream>
#include <stdio.h>
#include <map>
#include <vector>
#include <time.h>
#include <algorithm>
#include <math.h>

#include <pcl/console/print.h>
#include <pcl/console/parse.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>


#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/voxel_grid.h>
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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;
using namespace pcl;

pcl::PointCloud<pcl::PointXYZRGB>::Ptr  extract_small_set(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, Eigen::Vector3d point);
Eigen::Matrix4f estimate_contact_data(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, float *principle_curvature);


void show(pcl::PointCloud<pcl::PointXYZRGB>::Ptr showpc)
{

    pcl::visualization::PCLVisualizer vis3("VOXELIZED SAMPLES CLOUD");
    vis3.addPointCloud<pcl::PointXYZRGB>(showpc, "sample_cloud");

    vis3.setBackgroundColor(0, 0, 0); // 窗口背景色，默认[0,0,0]，范围[0~255,0~255,0~255]
    vis3.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample_cloud");			//设置点的大小，默认 1
    // vis3.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 1, 1, "sample_cloud");
    // vis3.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.5, "sample_cloud");			//设置点云透明度，默认 1 【Float going from 0.0 (transparent) to 1.0 (opaque)】
    vis3.addCoordinateSystem(100.0);
    // vis3.initCameraParameters ();
    vis3.spin();
}

void trans2center(pcl::PointCloud<pcl::PointXYZRGB> &cloud)
{

    Eigen::Matrix4f transform = Eigen::Matrix4f::Identity();
    Eigen::Matrix3f covariance_matrix;
    Eigen::Vector4f xyz_centroid_get;
    Eigen::Vector3f xyz_centroid, translation;

    pcl::compute3DCentroid(cloud, xyz_centroid_get);
    for (int i = 0; i < 3; i++)
    {
        xyz_centroid[i] = xyz_centroid_get[i];
    }
    pcl::computeCovarianceMatrix(cloud, xyz_centroid_get, covariance_matrix);

    Eigen::EigenSolver<Eigen::Matrix3f> eigensolver(covariance_matrix);
    Eigen::VectorXf eigen_values = eigensolver.eigenvalues().real();
    Eigen::Matrix3f eigen_vectors = eigensolver.eigenvectors().real();

    translation = -eigen_vectors.transpose() * xyz_centroid;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            transform(i, j) = (eigen_vectors.transpose())(i, j);
        }
    }
    for (int i = 0; i < 3; i++)
    {
        transform(i, 3) = translation[i];
    }
    pcl::transformPointCloud(cloud, cloud, transform);
}


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

    // Read in the cloud data
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::io::loadPCDFile<pcl::PointXYZRGB>(argv[pcd_file_indices[0]], *cloud);
    //Set color of point cloud
	size_t max_i = cloud->points.size();
	for (size_t i = 0; i < max_i; i++)
	{
		cloud->points[i].r = (uint8_t)255;
		cloud->points[i].g = (uint8_t)255;
		cloud->points[i].b = (uint8_t)255;
	}

    trans2center(*cloud);

    pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointNormal> mls_points;
    pcl::MovingLeastSquares<pcl::PointXYZRGB, pcl::PointXYZRGB> mls;
    // mls.setComputeNormals (true);
    // Set parameters
    mls.setInputCloud (cloud);
    mls.setPolynomialOrder (3);
    mls.setSearchMethod (tree);
    mls.setSearchRadius (5);
    // Reconstruct
    mls.setCacheMLSResults(true);
    mls.getCacheMLSResults();
    mls.process(*cloud1);
    cloud = cloud1;
    show(cloud);
   

    return 0;
}

//備份
/*
void path_generater::dynamic_adjust_path(Spline *origin_path, Spline pre_path)
{

    Spline *boundary = new Spline();
    if(compute_boundary(pre_path, boundary) == 0){
        printf("generate boundary fail\n");
        return;  
    } 
    
    std::map<double, std::vector<double>> new_path;
    double miny = origin_path->miny(), maxy = origin_path->bigy(), dy = 0;
    Eigen::Vector3f area_point, com_norm;
    Eigen::Vector3d path_node, bound_point;
    std::vector<float> Distance;
    pcl::Indices pointIdx;
    pcl::PointXYZRGB point;

    int NumOfNode = 10;  
    for (int i = 0; i < NumOfNode; i++){

        //注意數量有是否正確
        dy = ((maxy - miny) / NumOfNode * i) + miny;
        path_node = origin_path->point(dy);

        //use area_point to fix path node
        //area point is origin_path node's down bounded point
        area_point = Area2Cloud(path_node, 0, 1);
        
        /// check if overflow
        if (area_point[1] < boundary->miny() || area_point[1] > boundary->bigy()){
            point.x = path_node[0];
            point.y = path_node[1];
            point.z = path_node[2];
            kdtree.nearestKSearch(point, 3, pointIdx, Distance);
            std::vector<double> contain = {cloud->points[pointIdx[0]].x, cloud->points[pointIdx[0]].z};        
            new_path[cloud->points[pointIdx[0]].y] = contain;
            continue;
        }
            
        bound_point = boundary->point(area_point[1]);
   
        com_norm << (area_point[0] - bound_point[0]),
            (area_point[1] - bound_point[1]),
            (area_point[2] - bound_point[2]);

        if ((area_point[0] - bound_point[0]) > 0)
        {

            path_node << path_node[0] - com_norm[0],
                path_node[1] - com_norm[1],
                path_node[2] - com_norm[2];
        }
        else
        {
            path_node << path_node[0] + com_norm[0],
                path_node[1] + com_norm[1],
                path_node[2] + com_norm[2];
        }
        if(isnan(path_node[0])){
            std::cout << "fix path node NAN" << std::endl;
            continue;
        }
        point.x = path_node[0];
        point.y = path_node[1];
        point.z = path_node[2];
        kdtree.nearestKSearch(point, 5, pointIdx, Distance);

        std::vector<double> contain = {cloud->points[pointIdx[0]].x, cloud->points[pointIdx[0]].z};        
        new_path[cloud->points[pointIdx[0]].y] = contain;

    }
 
    int node_number = new_path.size(), index = 0;
    double point_x[node_number], point_y[node_number], point_z[node_number];
    for(auto& i: new_path){
        point_x[index] = i.second[0];
        point_y[index] = i.first;
        point_z[index] = i.second[1];
        index++;
        // std::cout << i.second[0] << " "<< i.first << " "<< i.second[1] << std::endl;
    }

    origin_path->restart(index, point_y, point_x, point_z);
    drawpath(*boundary,10,255,10);
    delete boundary;
    boundary = NULL;
}
*/