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
#include <pcl/features/principal_curvatures.h>
#include <pcl/filters/extract_indices.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;
using namespace pcl;

pcl::PointCloud<pcl::PointXYZRGB>::Ptr  extract_small_set(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, Eigen::Vector3d point);
Eigen::Matrix4f estimate_contact_data(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, float *principle_curvature);
Eigen::Matrix4f
compute_transform(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, pcl::PointXYZRGB point, float *principle_curvature);

pcl::KdTreeFLANN<pcl::PointXYZRGB> kdtree;
pcl::PointCloud<pcl::Normal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::Normal>);
pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principal_curvatures (new pcl::PointCloud<pcl::PrincipalCurvatures> ());
// Setup the principal curvatures computation


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

std::vector<int> 
rangedX_index(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, int position){
    std::vector<int> indices;
    pcl::PassThrough<pcl::PointXYZRGB> ptfilter (true);
    ptfilter.setInputCloud (cloud);
    ptfilter.setFilterFieldName ("x");
    ptfilter.setFilterLimits (-2+position, 2+position);
    ptfilter.filter(indices);

    return indices;
}

std::map<double, std::vector<double>> 
insert_point(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, 
    std::vector<int> indices, Eigen::Vector3f PlanePoint){ 

    std::vector<int> El, Er;
    Eigen::Vector3f normal(1,0,0);                                 
    
    //distinguish right or left points of plane
    for (int &i: indices) {
        Eigen::Vector3f point;
        point[0] = cloud->points[i].x;
        point[1] = cloud->points[i].y;
        point[2] = cloud->points[i].z;
        float distance2plane = (point-PlanePoint).dot(normal);
        if(distance2plane > 0){
            El.push_back(i);
        } 
        else if(distance2plane < 0){
            Er.push_back(i);
        }
	} 
    int El_flag[int(El.size())] = {0};
    int Er_flag[int(Er.size())] = {0};

    //get point pair
    Eigen::Vector3f Vector(0,0,0), point_l, point_r;
    std::map<float,int> compare;
    std::vector<int> left_pair, right_pair;
    int rp_index = 0;
    for(int i = 0 ; i < int(El.size()); i++){
        if(El_flag[i] == 0){ 
            //find the clostest point of right set betwen left_point
            compare.clear();
            for(int j = 0 ; j < int(Er.size()); j++){
                Vector[0] = cloud->points[El[i]].x - cloud->points[Er[j]].x;
                Vector[1] = cloud->points[El[i]].y - cloud->points[Er[j]].y;
                Vector[2] = cloud->points[El[i]].z - cloud->points[Er[j]].z;
                //index of Er at value of map 
                compare[Vector.norm()] = j;   
            }
            if(Er_flag[compare.begin()->second] == 0){
                //index of point cloud
                rp_index = Er[compare.begin()->second];
                right_pair.push_back(rp_index);
                //index of Er
                Er_flag[compare.begin()->second] = 1;
            }
            else 
                continue;

            //find the clostest point of left set betwen right_point which is clostest point
            compare.clear();
            for(int j = 0 ; j < int(El.size()); j++){
                Vector[0] = cloud->points[rp_index].x - cloud->points[El[j]].x;
                Vector[1] = cloud->points[rp_index].y - cloud->points[El[j]].y;
                Vector[2] = cloud->points[rp_index].z - cloud->points[El[j]].z;
                compare[Vector.norm()] = j;   
            }
            if(El_flag[compare.begin()->second] == 0){
                //index of point cloud 
                left_pair.push_back(El[compare.begin()->second]);
                //index of El
                El_flag[compare.begin()->second] = 1;
            }
        }
    }

    //insert new point to cloud
    pcl::PointCloud<pcl::PointXYZRGB> insert_cloud;
    insert_cloud.width = left_pair.size();
	insert_cloud.height = 1;
	insert_cloud.points.resize(insert_cloud.width*insert_cloud.height);
    //use for spline
    std::map<double, std::vector<double>> Node;
  
    for(int i=0; i < left_pair.size(); i++){
        int index_right = right_pair[i], index_left = left_pair[i];
        float t = (PlanePoint[0] - cloud->points[index_right].x)/(cloud->points[index_left].x - cloud->points[index_right].x);
        insert_cloud.points[i].x = PlanePoint[0];
		insert_cloud.points[i].y = cloud->points[index_right].y + t*(cloud->points[index_left].y - cloud->points[index_right].y);
		insert_cloud.points[i].z = cloud->points[index_right].z + t*(cloud->points[index_left].z - cloud->points[index_right].z);
        insert_cloud.points[i].r = (uint8_t)255;
        insert_cloud.points[i].g = (uint8_t)0;
        insert_cloud.points[i].b = (uint8_t)0;

        std::vector<double> temp = {insert_cloud.points[i].x,  insert_cloud.points[i].z};
        Node[insert_cloud.points[i].y] = temp;
    }
    // *cloud = (*cloud) + insert_cloud;
    return Node;
}

void estimate_normal(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud){

    pcl::NormalEstimation<pcl::PointXYZRGB, pcl::Normal> normal_estimation;
    normal_estimation.setInputCloud (cloud);
    pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZRGB>);
    normal_estimation.setSearchMethod (tree);
    // pcl::PointCloud<pcl::Normal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::Normal>);
    normal_estimation.setRadiusSearch(3);
    // normal_estimation.setKSearch (10); 
    normal_estimation.compute (*cloud_with_normals);

    // Setup the principal curvatures computation
    // pcl::PrincipalCurvaturesEstimation<pcl::PointXYZRGB, pcl::Normal, pcl::PrincipalCurvatures> principal_curvatures_estimation;

    // principal_curvatures_estimation.setInputCloud (cloud);
    // principal_curvatures_estimation.setInputNormals (cloud_with_normals);
    // principal_curvatures_estimation.setSearchMethod (tree);
    // // principal_curvatures_estimation.setKSearch (10);  
    // principal_curvatures_estimation.setRadiusSearch(2);
    // principal_curvatures_estimation.compute (*principal_curvatures);
}
/*
pcl::PointCloud<pcl::PointXYZRGB>::Ptr 
extract_small_set(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, Eigen::Vector3d point){
    
    //Set the point 
    pcl::PointXYZRGB searchPoint;
    searchPoint.x = point[0];
    searchPoint.y = point[1];
    searchPoint.z = point[2];
    pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZRGB>);
    //Set kd-tree for point set
    tree->setInputCloud(cloud);
    float radius = 100;
    // std::vector<int> indices;
    pcl::Indices indices;
    std::vector<float> distance;
    tree->radiusSearch(searchPoint, radius, indices, distance);

    pcl::IndicesPtr  index_ptr(new std::vector<int>(indices)); 
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr extra_cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::ExtractIndices<pcl::PointXYZRGB> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(index_ptr);
    extract.setNegative(false);//如果设为true,可以提取指定index之外的点云
    extract.filter(*extra_cloud);
    extra_cloud->points[0] = searchPoint;

    return extra_cloud;
}
*/
//Small set of point cloud use to estimate contact data 
/*
Eigen::Matrix4f estimate_contact_data(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, float *principle_curvature){
    pcl::NormalEstimation<pcl::PointXYZRGB, pcl::Normal> normal_estimation;
    normal_estimation.setInputCloud (cloud);
    pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZRGB>);
    normal_estimation.setSearchMethod (tree);
    pcl::PointCloud<pcl::Normal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::Normal>);
    normal_estimation.setRadiusSearch(1);
    // normal_estimation.setKSearch (10); 
    normal_estimation.compute (*cloud_with_normals);


    // Setup the principal curvatures computation
    pcl::PrincipalCurvaturesEstimation<pcl::PointXYZRGB, pcl::Normal, pcl::PrincipalCurvatures> principal_curvatures_estimation;

    principal_curvatures_estimation.setInputCloud (cloud);
    principal_curvatures_estimation.setInputNormals (cloud_with_normals);
    principal_curvatures_estimation.setSearchMethod (tree);
    principal_curvatures_estimation.setKSearch (10);  
    // principal_curvatures_estimation.setRadiusSearch(1);

    pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principal_curvatures (new pcl::PointCloud<pcl::PrincipalCurvatures> ());
    principal_curvatures_estimation.compute (*principal_curvatures);
    principle_curvature[0] = principal_curvatures->points[0].pc1;
    principle_curvature[1] = principal_curvatures->points[0].pc2;

    Eigen::Vector3f curvatureVector, normalVector;

    curvatureVector <<  principal_curvatures->points[0].principal_curvature_x,
                        principal_curvatures->points[0].principal_curvature_y,
                        principal_curvatures->points[0].principal_curvature_z;

    normalVector << cloud_with_normals->points[0].normal_x,
                    cloud_with_normals->points[0].normal_y,
                    cloud_with_normals->points[0].normal_z; 
    // printf("normal \n(%f, %f,%f)\n",normalVector[0], normalVector[1], normalVector[2]);
    Eigen::Vector3f crossVector = normalVector.cross(curvatureVector);
    // printf("Small eigein values's vector \n(%f, %f,%f)\n",crossVector[0],crossVector[1],crossVector[2]);

    // noa Rotation Matrix
    //bigger curvature dir, smallerr curvature dir, normal, point
    Eigen::Matrix4f Y; Y << crossVector(0), curvatureVector(0), normalVector(0), cloud->points[0].x,
                            crossVector(1), curvatureVector(1), normalVector(1), cloud->points[0].y,
                            crossVector(2), curvatureVector(2), normalVector(2), cloud->points[0].z,
                                            0 ,             0 ,              0 ,                   1 ;
    Eigen::Matrix4f X; X << Eigen::Matrix<float, 4, 4>::Identity();
    Eigen::Matrix4f pt2Base = X.transpose() * Y;

    return pt2Base;
}
*/

/*
void AddArea2Cloud(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, Eigen::Vector3d point){

    double toolRadius = 20, depth = 0.1;
    Eigen::Matrix4f transform;
    float principle_curvature[2];
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr extra_cloud;
    extra_cloud = extract_small_set(cloud, point);
    transform = estimate_contact_data(extra_cloud, principle_curvature);

    // printf("curvature1 : %lf\n",principle_curvature[0]);
    // printf("curvature2 : %lf\n",principle_curvature[1]);
    
    double longAxis, shortAxis;
    longAxis = std::sqrt(std::pow(1/principle_curvature[1],2)-std::pow(std::abs(1/principle_curvature[1])-depth,2));
    if(longAxis>toolRadius){
        longAxis =toolRadius;
    }
    shortAxis = std::sqrt(std::pow(1/principle_curvature[0],2)-std::pow(std::abs(1/principle_curvature[0])-depth,2));
    if(shortAxis>toolRadius){
        shortAxis =toolRadius;
    }
    // printf("longAxis : %lf\n",longAxis);
    // printf("shortAxis : %lf\n",shortAxis);

    //draw a Ellipse in origin ,and trasform to path node 
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr ellipse (new pcl::PointCloud<pcl::PointXYZRGB>);
    std::uint8_t r(255), g(0), b(0);
    for (float angle(0.0); angle <= 360.0; angle += 1){
        pcl::PointXYZRGB basic_point;
        basic_point.x = longAxis * std::cos (pcl::deg2rad(angle));
        basic_point.y = shortAxis * std::sin (pcl::deg2rad(angle));
        basic_point.z = 0;
        basic_point.r = 255;
        basic_point.g = 0;
        basic_point.b = 0;
        ellipse->points.push_back(basic_point);
    }
    pcl::transformPointCloud(*ellipse, *ellipse, transform);
    *cloud += *ellipse;

    
    // std::cout << curvatureVector << std::endl<< std::endl;
    // std::cout << normalVector << std::endl<< std::endl;
    // std::cout << pt2Base << std::endl<< std::endl;
}
*/

void Area2Cloud(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, Eigen::Vector3d point){

    double toolRadius = 20, depth = 0.1;
    Eigen::Matrix4f transform;
    float principle_curvature[2];

    // kdtree.setInputCloud (cloud);
    pcl::PointXYZRGB SearchPoint;
    SearchPoint.x = point[0];
    SearchPoint.y = point[1];
    SearchPoint.z = point[2];

    transform = compute_transform(cloud, SearchPoint, principle_curvature);
    
    double longAxis, shortAxis;
    longAxis = std::sqrt(std::pow(1/principle_curvature[1],2)-std::pow(std::abs(1/principle_curvature[1])-depth,2));
    if(longAxis>toolRadius){
        longAxis =toolRadius;
    }
    shortAxis = std::sqrt(std::pow(1/principle_curvature[0],2)-std::pow(std::abs(1/principle_curvature[0])-depth,2));
    if(shortAxis>toolRadius){
        shortAxis =toolRadius;
    }

    //draw a Ellipse in origin ,and trasform to path node 
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr ellipse (new pcl::PointCloud<pcl::PointXYZRGB>);
    std::uint8_t r(255), g(0), b(0);
    for (float angle(0.0); angle <= 360.0; angle += 1){
        pcl::PointXYZRGB basic_point;
        basic_point.x = longAxis * std::cos (pcl::deg2rad(angle));
        basic_point.y = shortAxis * std::sin (pcl::deg2rad(angle));
        basic_point.z = 0;
        basic_point.r = 255;
        basic_point.g = 0;
        basic_point.b = 0;
        ellipse->points.push_back(basic_point);
    }
    pcl::transformPointCloud(*ellipse, *ellipse, transform);
    *cloud += *ellipse;

    
    // std::cout << curvatureVector << std::endl<< std::endl;
    // std::cout << normalVector << std::endl<< std::endl;
    // std::cout << pt2Base << std::endl<< std::endl;
}
/*
void spline_path(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, Eigen::Vector3f plane_point){
    struct timeval mstart, mend;
    unsigned  long times;

    gettimeofday(&mstart,NULL);
    Eigen::Vector3f normal(1,0,0);
    // Eigen::Vector3f plane_point(90,0,0);
    std::map<double, std::vector<double>> Node;
    Node = insert_point(cloud, rangedX_index(cloud, plane_point[0]), plane_point);
    gettimeofday(&mend,NULL);
    times = 1000000 * (mend.tv_sec-mstart.tv_sec)+ mend.tv_usec-mstart.tv_usec;
    printf("insert_point time: %lu \n",times);

    int node_number = Node.size(), index = 0;
    double point_x[node_number], point_y[node_number], point_z[node_number];
    for(auto& i: Node){
        point_x[index] = i.second[0];
        point_y[index] = i.first;
        point_z[index] = i.second[1];
        index++;
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    // const gsl_interp_type *Method = gsl_interp_cspline_periodic;
    const gsl_interp_type *Method = gsl_interp_steffen;
    gsl_spline *splineYX = gsl_spline_alloc (Method, node_number);
    gsl_spline *splineYZ = gsl_spline_alloc (Method, node_number);

    gsl_spline_init (splineYX, point_y, point_x, node_number);
    gsl_spline_init (splineYZ, point_y, point_z, node_number);
    double dy;
    Eigen::Vector3d point;
    // Eigen::Vector3f noraml;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr extra_cloud;

    // point << point_x[5], point_y[5], point_z[5];
    // AddArea2Cloud(cloud, point);
    // float principle_curvature[2];
    // int increase = index/20;


    for (int i = 0; i < 9; i++){
        gettimeofday(&mstart,NULL);
        // dy = (1 - i / 10) * point_y[0] + (i / 10) * point_y[node_number-1];
        dy = ((point_y[node_number-1]-point_y[0])/10*i)+point_y[0];
        // printf("dy is %lf \n",point_y[0]-point_y[node_number-1]);
        // printf("dy is %lf \n",((point_y[node_number-1]-point_y[0])/10*i)+point_y[0]);
        // printf("dy is %lf \n",dy);
        // printf("dy is %lf \n",point_y[0]);
        // std::cout << "x:" << gsl_spline_eval (splineYX, dy, acc) << "  ";
        
        // std::cout << "ori_y:" << point_y[i] << "  ";
        // std::cout << "ori_z:" << point_z[i] << "  ";
        // dy = (point_y[i]-point_y[0])/(point_y[node_number-1]-point_y[0]);
        // cout << dy << " ";
        // std::cout << "z1:" << gsl_spline_eval (splineYZ, point_y[i]+0.02, acc) << "  ";
        // std::cout << "z2:" << gsl_spline_eval (splineYZ, dy, acc) << "  ";
        
        // std::cout << std::endl;
        point << gsl_spline_eval(splineYX, dy, acc), dy, gsl_spline_eval(splineYZ, dy, acc);
        // point << point_x[i], point_y[i], point_z[i];
        // printf("Points(%f, %f, %f)\n",gsl_spline_eval(splineYX, dy, acc), dy, gsl_spline_eval(splineYZ, dy, acc));
        extra_cloud = extract_small_set(cloud, point);
        // printf("point(%f, %f, %f)\n",point[0],point[1],point[2]);
        // printf("norma(%f, %f, %f)\n",normal[0],normal[1],normal[2]);
        // estimate_contact_data(extra_cloud, principle_curvature);

        AddArea2Cloud(cloud, point);
        // printf("\n");
        gettimeofday(&mend,NULL);
        times = 1000000 * (mend.tv_sec-mstart.tv_sec)+ mend.tv_usec-mstart.tv_usec;
        printf("estimat contact use time: %lu \n",times);
    }
    printf("Get a path !\n");
    gsl_spline_free (splineYX);
    gsl_spline_free (splineYZ);
    gsl_interp_accel_free (acc);
}
*/


Eigen::Matrix4f
compute_transform(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, pcl::PointXYZRGB point, float *principle_curvature){
    Eigen::Vector3f curvatureVector, normalVector;
    // principal_curvatures_estimation.setInputCloud (cloud);
    // principal_curvatures_estimation.setInputNormals (cloud_with_normals);
    // principal_curvatures_estimation.setSearchMethod (tree);
    // principal_curvatures_estimation.setKSearch (10);  
    // principal_curvatures_estimation.setRadiusSearch(2);
    // pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principal_curvatures (new pcl::PointCloud<pcl::PrincipalCurvatures> ());
    // principal_curvatures_estimation.compute (*principal_curvatures);

    // int K = 1;                             //搜索最近邻的点数
    // kdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance)
    // std::vector<int> pointIdx;
    

    
    std::vector<float> Distance;
    pcl::Indices pointIdx;

    float radius = 1;
    kdtree.radiusSearch (point, radius, pointIdx, Distance);

    pcl::PrincipalCurvaturesEstimation<pcl::PointXYZRGB, pcl::Normal, pcl::PrincipalCurvatures> principal_curvatures_estimation;
    // float vector[3];
    principal_curvatures_estimation.computePointPrincipalCurvatures(*cloud_with_normals, 
                                                                    pointIdx[0], 
                                                                    pointIdx, 
                                                                    curvatureVector[0], 
                                                                    curvatureVector[1], 
                                                                    curvatureVector[2],
                                                                    principle_curvature[0],
                                                                    principle_curvature[1]
                                                                    );
       
    // principle_curvature[0] = principal_curvatures->points[pointIdx[0]].pc1;
    // principle_curvature[1] = principal_curvatures->points[pointIdx[0]].pc2;
    // curvatureVector <<  vector[0], vector[1], vector[2];
     
    // curvatureVector <<  principal_curvatures->points[pointIdx[0]].principal_curvature_x,
    //                     principal_curvatures->points[pointIdx[0]].principal_curvature_y,
    //                     principal_curvatures->points[pointIdx[0]].principal_curvature_z;

    // std::cout << cloud_with_normals->points[0].normal_x << std::endl;
    normalVector << cloud_with_normals->points[pointIdx[0]].normal_x,
                    cloud_with_normals->points[pointIdx[0]].normal_y,
                    cloud_with_normals->points[pointIdx[0]].normal_z; 

    // printf("normal \n(%f, %f,%f)\n",normalVector[0], normalVector[1], normalVector[2]);
    Eigen::Vector3f crossVector = normalVector.cross(curvatureVector);
    // printf("Small eigein values's vector \n(%f, %f,%f)\n",crossVector[0],crossVector[1],crossVector[2]);

    // noa Rotation Matrix
    //bigger curvature dir, smallerr curvature dir, normal, point
    Eigen::Matrix4f Y; Y << crossVector(0), curvatureVector(0), normalVector(0), point.x,
                            crossVector(1), curvatureVector(1), normalVector(1), point.y,
                            crossVector(2), curvatureVector(2), normalVector(2), point.z,
                                            0 ,             0 ,              0 ,                   1 ;
    Eigen::Matrix4f X; X << Eigen::Matrix<float, 4, 4>::Identity();
    Eigen::Matrix4f pt2Base = X.transpose() * Y;

    return pt2Base;

}


void path_track(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, Eigen::Vector3f plane_point){

    struct timeval mstart, mend;
    unsigned  long times;
    gettimeofday(&mstart,NULL);
    //Insert Point Node with a line
    Eigen::Vector3f normal(1,0,0);
    std::map<double, std::vector<double>> Node;
    Node = insert_point(cloud, rangedX_index(cloud, plane_point[0]), plane_point);
    gettimeofday(&mend,NULL);
    times = 1000000 * (mend.tv_sec-mstart.tv_sec)+ mend.tv_usec-mstart.tv_usec;
    printf("insert_point time: %lu \n",times);

    int node_number = Node.size(), index = 0;
    double point_x[node_number], point_y[node_number], point_z[node_number];
    for(auto& i: Node){
        point_x[index] = i.second[0];
        point_y[index] = i.first;
        point_z[index] = i.second[1];
        index++;
    }
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    // const gsl_interp_type *Method = gsl_interp_cspline_periodic;
    const gsl_interp_type *Method = gsl_interp_steffen;
    gsl_spline *splineYX = gsl_spline_alloc (Method, node_number);
    gsl_spline *splineYZ = gsl_spline_alloc (Method, node_number);

    gsl_spline_init (splineYX, point_y, point_x, node_number);
    gsl_spline_init (splineYZ, point_y, point_z, node_number);
    double dy;
    Eigen::Vector3d point;
    // Eigen::Vector3f noraml;

    for (int i = 0; i < 20; i++){
        gettimeofday(&mstart,NULL);

        dy = ((point_y[node_number-1]-point_y[0])/20*i)+point_y[0];
        point << gsl_spline_eval(splineYX, dy, acc), dy, gsl_spline_eval(splineYZ, dy, acc);
        Area2Cloud(cloud, point);
        // printf("\n");
        gettimeofday(&mend,NULL);
        times = 1000000 * (mend.tv_sec-mstart.tv_sec)+ mend.tv_usec-mstart.tv_usec;
        printf("estimat contact use time: %lu \n",times);
    }
    printf("Get a path !\n");
    gsl_spline_free (splineYX);
    gsl_spline_free (splineYZ);
    gsl_interp_accel_free (acc);
}

void slicing_method(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud)
{
    printf("Start Path Planning!\n");

    struct timeval start, end;
    unsigned long time;
    gettimeofday(&start, NULL);
    kdtree.setInputCloud (cloud);
    // Plane model
    Eigen::Vector3f normal(1, 0, 0);
    Eigen::Vector3f plane_point(100, 0, 0);
    std::vector<int> indices;
    // find the bounding point
    pcl::PointXYZRGB min_pt, max_pt;
    pcl::getMinMax3D(*cloud, min_pt, max_pt);


    int step_size = 40, path_nums = 0;
    min_pt.x += step_size;
    while (min_pt.x < max_pt.x)
    {

        plane_point[0] = min_pt.x;
        path_track(cloud, plane_point);
        // spline_path(cloud, plane_point);
        min_pt.x += step_size;
        path_nums++;

        // cout << "use time: " << time << endl;
    }
    // spline_path(cloud, plane_point);
    gettimeofday(&end, NULL);
    time = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    cout << "use time: " << time << endl;
    cout << "number of paths: " << path_nums << endl;
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
    pcl::io::loadPCDFile<pcl::PointXYZRGB>(argv[pcd_file_indices[0]], *cloud);
    // Create the filtering object
    // pcl::VoxelGrid<pcl::PointXYZRGB> sor;
    // sor.setInputCloud(cloud);
    // sor.setLeafSize(0.1, 0.1, 0.1);
    // sor.filter(*cloud);
    //Set color of point cloud
	size_t max_i = cloud->points.size();
	for (size_t i = 0; i < max_i; i++)
	{
		cloud->points[i].r = (uint8_t)255;
		cloud->points[i].g = (uint8_t)255;
		cloud->points[i].b = (uint8_t)255;
	}

    trans2center(*cloud);
    // show(cloud);
    // slicing_method(cloud);
    estimate_normal(cloud);
    slicing_method(cloud);
    show(cloud);
   

    return 0;
}
