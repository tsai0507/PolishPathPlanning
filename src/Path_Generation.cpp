#include "Path_Generate.h"

path_generater::path_generater(std::string cloud_name, double Radius) : toolRadius(Radius), file_name(cloud_name)
{
    cloud = pcl::PointCloud<pcl::PointXYZRGB>::Ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
    other_cloud = pcl::PointCloud<pcl::PointXYZRGB>::Ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
    cloud_with_normals = pcl::PointCloud<pcl::Normal>::Ptr(new pcl::PointCloud<pcl::Normal>);
    if (pcl::io::loadPCDFile<pcl::PointXYZRGB>(cloud_name, *cloud) == -1)
    {
        PCL_ERROR("Cloudn't read file!\n");
    }
    else
    {
        Path_set = {};
        size_t max_i = cloud->points.size();

        // pcl::StatisticalOutlierRemoval<pcl::PointXYZRGB> sor;
        // sor.setInputCloud (cloud);
        // sor.setMeanK (50);
        // sor.setStddevMulThresh (1.0);
        // sor.filter (*cloud);

        for (size_t i = 0; i < max_i; i++)
        {
            cloud->points[i].r = (uint8_t)255;
            cloud->points[i].g = (uint8_t)255;
            cloud->points[i].b = (uint8_t)255;
            // cloud->points[i].x *= 1000;
            // cloud->points[i].y *= 1000;
            // cloud->points[i].z *= 1000;
        }
        transform = Eigen::Matrix4f::Identity();
        normal << 1, 0, 0;
    }
};

void path_generater::show()
{
    *other_cloud += *cloud;
    pcl::visualization::PCLVisualizer viewer("VOXELIZED SAMPLES CLOUD");
    viewer.addPointCloud<pcl::PointXYZRGB>(other_cloud, "sample_cloud");

    viewer.setBackgroundColor(0, 0, 0);                                                                        // 窗口背景色，默认[0,0,0]，范围[0~255,0~255,0~255]
    viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample_cloud"); // 设置点的大小，默认 1
    // viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 1, 1, "sample_cloud");
    // viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, 0.5, "sample_cloud");			//设置点云透明度，默认 1 【Float going from 0.0 (transparent) to 1.0 (opaque)】
    viewer.addCoordinateSystem(100.0);
    // viewer.initCameraParameters ();

    viewer.spin();
}

void path_generater::voxel_down(const float x, const float y, const float z)
{
    pcl::VoxelGrid<pcl::PointXYZRGB> sor;
    sor.setInputCloud(cloud);
    sor.setLeafSize(x, y, z);
    sor.filter(*cloud);
}

void path_generater::trans2center()
{

    Eigen::Matrix3f covariance_matrix;
    Eigen::Vector4f xyz_centroid_get;
    Eigen::Vector3f xyz_centroid, translation;

    pcl::compute3DCentroid(*cloud, xyz_centroid_get);
    for (int i = 0; i < 3; i++)
    {
        xyz_centroid[i] = xyz_centroid_get[i];
    }
    pcl::computeCovarianceMatrix(*cloud, xyz_centroid_get, covariance_matrix);

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
    pcl::transformPointCloud(*cloud, *cloud, transform);
}

std::vector<int> path_generater::rangedX_index(int position)
{
    std::vector<int> indices;
    pcl::PassThrough<pcl::PointXYZRGB> ptfilter(true);
    ptfilter.setInputCloud(cloud);
    ptfilter.setFilterFieldName("x");
    ptfilter.setFilterLimits(-2 + position, 2 + position);
    ptfilter.filter(indices);

    return indices;
}


std::map<double, std::vector<double>>
path_generater::insert_point(std::vector<int> indices, Eigen::Vector3f PlanePoint)
{
    std::vector<int> El, Er;

    // distinguish right or left points of plane
    for (int &i : indices)
    {
        Eigen::Vector3f point;
        point[0] = cloud->points[i].x;
        point[1] = cloud->points[i].y;
        point[2] = cloud->points[i].z;
        float distance2plane = (point - PlanePoint).dot(normal);
        if (distance2plane > 0)
        {
            El.push_back(i);
        }
        else if (distance2plane < 0)
        {
            Er.push_back(i);
        }
    }
    int El_flag[int(El.size())] = {0};
    int Er_flag[int(Er.size())] = {0};

    // get point pair
    Eigen::Vector3f Vector(0, 0, 0), point_l, point_r;
    std::map<float, int> compare;
    std::vector<int> left_pair, right_pair;
    int rp_index = 0;
    for (int i = 0; i < int(El.size()); i++)
    {
        if (El_flag[i] == 0)
        {
            // find the clostest point of right set betwen left_point
            compare.clear();
            for (int j = 0; j < int(Er.size()); j++)
            {
                Vector[0] = cloud->points[El[i]].x - cloud->points[Er[j]].x;
                Vector[1] = cloud->points[El[i]].y - cloud->points[Er[j]].y;
                Vector[2] = cloud->points[El[i]].z - cloud->points[Er[j]].z;
                // index of Er at value of map
                compare[Vector.norm()] = j;
            }
            if (Er_flag[compare.begin()->second] == 0)
            {
                // index of point cloud
                rp_index = Er[compare.begin()->second];
                right_pair.push_back(rp_index);
                // index of Er
                Er_flag[compare.begin()->second] = 1;
            }
            else
                continue;

            // find the clostest point of left set betwen right_point which is clostest point
            compare.clear();
            for (int j = 0; j < int(El.size()); j++)
            {
                Vector[0] = cloud->points[rp_index].x - cloud->points[El[j]].x;
                Vector[1] = cloud->points[rp_index].y - cloud->points[El[j]].y;
                Vector[2] = cloud->points[rp_index].z - cloud->points[El[j]].z;
                compare[Vector.norm()] = j;
            }
            if (El_flag[compare.begin()->second] == 0)
            {
                // index of point cloud
                left_pair.push_back(El[compare.begin()->second]);
                // index of El
                El_flag[compare.begin()->second] = 1;
            }
        }
    }

    // insert new point to cloud
    pcl::PointCloud<pcl::PointXYZRGB> insert_cloud;
    insert_cloud.width = left_pair.size();
    insert_cloud.height = 1;
    insert_cloud.points.resize(insert_cloud.width * insert_cloud.height);
    // use for spline
    std::map<double, std::vector<double>> Node;

    for (int i = 0; i < left_pair.size(); i++)
    {
        int index_right = right_pair[i], index_left = left_pair[i];
        float t = (PlanePoint[0] - cloud->points[index_right].x) / (cloud->points[index_left].x - cloud->points[index_right].x);
        insert_cloud.points[i].x = PlanePoint[0];
        insert_cloud.points[i].y = cloud->points[index_right].y + t * (cloud->points[index_left].y - cloud->points[index_right].y);
        insert_cloud.points[i].z = cloud->points[index_right].z + t * (cloud->points[index_left].z - cloud->points[index_right].z);
        insert_cloud.points[i].r = (uint8_t)0;
        insert_cloud.points[i].g = (uint8_t)0;
        insert_cloud.points[i].b = (uint8_t)255;

        std::vector<double> temp = {insert_cloud.points[i].x, insert_cloud.points[i].z};
        Node[insert_cloud.points[i].y] = temp;
    }
    // *cloud = (*cloud) + insert_cloud;
    *other_cloud = (*other_cloud) + insert_cloud;
    return Node;
}

// std::map<double, std::vector<double>>
// path_generater::insert_point(std::vector<int> indices, Eigen::Vector3f PlanePoint)
// {
//     std::vector<int> El, Er;
//     std::vector<int> left_pair, right_pair;
//     std::vector<int> idr(1), idl(1), idsave(1);
//     std::vector<float> dis(1);
    
//     // distinguish right or left points of plane
//     Eigen::Vector3f point;
//     for (int &i : indices)
//     {
//         point[0] = cloud->points[i].x;
//         point[1] = cloud->points[i].y;
//         point[2] = cloud->points[i].z;
//         float distance2plane = (point - PlanePoint).dot(normal);
//         if (distance2plane > 0)
//             El.push_back(i);
//         else if (distance2plane < 0)
//             Er.push_back(i);
//     }

//     PointCloudType::Ptr cloudEl(new PointCloudType);
//     PointCloudType::Ptr cloudEr(new PointCloudType);
//     pcl::KdTreeFLANN<PointType> treeEl;
//     pcl::KdTreeFLANN<PointType> treeEr;
//     pcl::copyPointCloud(*cloud, El, *cloudEl);
//     pcl::copyPointCloud(*cloud, Er, *cloudEr);
//     treeEl.setInputCloud(cloudEl);
//     treeEr.setInputCloud(cloudEr);

//     //left set which closet to right, right set which closet to left
//     for (int i = 0; i < int(El.size()); i++)
//     {
//         auto pl = cloudEl->points[i];
//         treeEr.nearestKSearch(pl, 1, idr, dis);

//         auto pr = cloudEr->points[idr[0]];
//         kdtree.nearestKSearch(pr, 1, idsave, dis);
//         right_pair.push_back(idsave[0]);
        
//         treeEl.nearestKSearch(pr, 1, idl, dis);
//         pl = cloudEl->points[idl[0]];
//         kdtree.nearestKSearch(pl, 1, idsave, dis);
//         left_pair.push_back(idsave[0]);
//     }

//     // insert new point to cloud
//     pcl::PointCloud<pcl::PointXYZRGB> insert_cloud;
//     insert_cloud.width = left_pair.size();
//     insert_cloud.height = 1;
//     insert_cloud.points.resize(insert_cloud.width * insert_cloud.height);
//     // use for spline
//     std::map<double, std::vector<double>> Node;

//     for (int i = 0; i < left_pair.size(); i++)
//     {
//         int index_right = right_pair[i], index_left = left_pair[i];
//         float t = (PlanePoint[0] - cloud->points[index_right].x) / (cloud->points[index_left].x - cloud->points[index_right].x);
//         insert_cloud.points[i].x = PlanePoint[0];
//         insert_cloud.points[i].y = cloud->points[index_right].y + t * (cloud->points[index_left].y - cloud->points[index_right].y);
//         insert_cloud.points[i].z = cloud->points[index_right].z + t * (cloud->points[index_left].z - cloud->points[index_right].z);
//         insert_cloud.points[i].r = (uint8_t)0;
//         insert_cloud.points[i].g = (uint8_t)0;
//         insert_cloud.points[i].b = (uint8_t)255;

//         std::vector<double> temp = {insert_cloud.points[i].x, insert_cloud.points[i].z};
//         Node[insert_cloud.points[i].y] = temp;
//     }
//     // *cloud = (*cloud) + insert_cloud;
//     *other_cloud = (*other_cloud) + insert_cloud;
//     return Node;
// }

void path_generater::slicing_method()
{
    struct timeval start, end;
    unsigned long time;
    gettimeofday(&start, NULL);

    // Plane model
    Eigen::Vector3f plane_point(0, 0, 0);
    std::vector<int> indices;
    // find the bounding point
    pcl::PointXYZRGB min_pt, max_pt;
    pcl::getMinMax3D(*cloud, min_pt, max_pt);

    int path_nums = 0, step_size = toolRadius * 2;
    min_pt.x += step_size/2;
    while (min_pt.x < max_pt.x)
    {
        plane_point[0] = min_pt.x;
        indices = rangedX_index(plane_point[0]);
        insert_point(indices, plane_point);
        min_pt.x += step_size;
        path_nums++;
    }

    gettimeofday(&end, NULL);
    time = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    std::cout << "use time: " << time << std::endl;
    std::cout << "number of paths: " << path_nums << std::endl;
    // pcl::io::savePCDFileASCII ("insert.pcd", *cloud);
}

void path_generater::estimate_normal()
{
    pcl::NormalEstimation<pcl::PointXYZRGB, pcl::Normal> normal_estimation;
    normal_estimation.setInputCloud(cloud);
    pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZRGB>);
    normal_estimation.setSearchMethod(tree);
    normal_estimation.setRadiusSearch(2.5);
    // normal_estimation.setKSearch (10);
    normal_estimation.compute(*cloud_with_normals);

}

void path_generater::Set_kdtree()
{
    kdtree.setInputCloud(cloud);
}

void path_generater::smooth(){

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointNormal> mls_points;
    pcl::MovingLeastSquares<pcl::PointXYZRGB, pcl::PointXYZRGB> mls;
    // mls.setComputeNormals (true);
    // Set parameters
    mls.setInputCloud (cloud);
    mls.setPolynomialOrder (3);
    mls.setSearchMethod (tree);
    mls.setSearchRadius (15);
    // Reconstruct
    // mls.setCacheMLSResults(true);
    // mls.getCacheMLSResults();
    mls.process(*cloud1);
    cloud = cloud1;
    std::string smooth_nam = "smooth_";
    smooth_nam += file_name;
    pcl::io::savePCDFileASCII(smooth_nam, *cloud1);
}

Eigen::Matrix4f
path_generater::compute_transform(pcl::PointXYZRGB point, float *principle_curvature)
{

    Eigen::Vector3f curvatureVector, normalVector;
    std::vector<float> Distance;
    pcl::Indices pointIdx;
    // float radius = 1;
    // kdtree.nearestKSearch(point, toolRadius*0.8, pointIdx, Distance);
    // kdtree.radiusSearch(point, toolRadius*0.8, pointIdx, Distance);
    kdtree.nearestKSearch(point, 10, pointIdx, Distance);

    pcl::PrincipalCurvaturesEstimation<pcl::PointXYZRGB, pcl::Normal, pcl::PrincipalCurvatures> principal_curvatures_estimation;
    principal_curvatures_estimation.computePointPrincipalCurvatures(*cloud_with_normals,
                                                                    pointIdx[0],
                                                                    pointIdx,
                                                                    curvatureVector[0],
                                                                    curvatureVector[1],
                                                                    curvatureVector[2],
                                                                    principle_curvature[0],
                                                                    principle_curvature[1]);
    normalVector << cloud_with_normals->points[pointIdx[0]].normal_x,
        cloud_with_normals->points[pointIdx[0]].normal_y,
        cloud_with_normals->points[pointIdx[0]].normal_z;
    Eigen::Vector3f crossVector = normalVector.cross(curvatureVector);

    // noa Rotation Matrix
    // bigger curvature dir, smallerr curvature dir, normal, point
    Eigen::Matrix4f Y;
    Y << crossVector(0), curvatureVector(0), normalVector(0), point.x,
        crossVector(1), curvatureVector(1), normalVector(1), point.y,
        crossVector(2), curvatureVector(2), normalVector(2), point.z,
        0, 0, 0, 1;
    Eigen::Matrix4f X;
    X << Eigen::Matrix<float, 4, 4>::Identity();
    Eigen::Matrix4f pt2Base = X.transpose() * Y;

    return pt2Base;
}

// flag draw the data, key up/down
Eigen::Vector3f
path_generater::Area2Cloud(Eigen::Vector3d point, bool flag, bool key)
{

    Eigen::Matrix4f transform;
    float principle_curvature[2];
    pcl::PointXYZRGB SearchPoint;
    SearchPoint.x = point[0];
    SearchPoint.y = point[1];
    SearchPoint.z = point[2];
    transform = compute_transform(SearchPoint, principle_curvature);
    double longAxis, shortAxis;
    if((principle_curvature[0] >= 0) && (principle_curvature[1] >= 0)){
        longAxis = std::sqrt(std::pow(1 / principle_curvature[1], 2) - std::pow(std::abs(1 / principle_curvature[1]) - depth, 2));
        if (longAxis > toolRadius){
            longAxis = toolRadius;
        }

        shortAxis = std::sqrt(std::pow(1 / principle_curvature[0], 2) - std::pow(std::abs(1 / principle_curvature[0]) - depth, 2));
        if (shortAxis > toolRadius){
            shortAxis = toolRadius;
        }
    }
    else{
        longAxis = std::abs(1 / principle_curvature[1]) - std::sqrt(std::pow(1 / principle_curvature[1], 2) - std::pow(toolRadius, 2));
        if (longAxis > toolthickness){
            longAxis = toolthickness;
        }

        shortAxis = std::abs(1 / principle_curvature[0]) - std::sqrt(std::pow(1 / principle_curvature[0], 2) - std::pow(toolRadius, 2));
        if (shortAxis > toolthickness){
            shortAxis = toolthickness;
        }
    }
    

    // draw a Ellipse in origin ,and trasform to path node
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr ellipse(new pcl::PointCloud<pcl::PointXYZRGB>);
    std::uint8_t r(255), g(0), b(0);
    for (float angle(0.0); angle <= 360.0; angle += 0.5)
    {
        pcl::PointXYZRGB basic_point;
        basic_point.x = longAxis * std::cos(pcl::deg2rad(angle));
        basic_point.y = shortAxis * std::sin(pcl::deg2rad(angle));
        basic_point.z = 0;
        basic_point.r = 255;
        basic_point.g = 0;
        basic_point.b = 0;
        ellipse->points.push_back(basic_point);
    }
    pcl::transformPointCloud(*ellipse, *ellipse, transform);

    // find the boundary point
    Eigen::Vector3f boundpoint;
    if (key == 0)
    {
        auto boundpoint_it = std::max_element(ellipse->begin(), ellipse->end(), [](auto &p1, auto &p2)
                                              { return p1.x < p2.x; });
        boundpoint << boundpoint_it->x, boundpoint_it->y, boundpoint_it->z;

        //use for compute coverage
        auto for_min = std::min_element(ellipse->begin(), ellipse->end(), [](auto &p1, auto &p2)
                                              { return p1.x < p2.x; });
        double comput_lan =   (for_min-> x - boundpoint_it->x)/2; 
        compute_coverage(point, comput_lan);                           
    }
    else
    {
        auto boundpoint_it = std::min_element(ellipse->begin(), ellipse->end(), [](auto &p1, auto &p2)
                                              { return p1.x < p2.x; });
        boundpoint << boundpoint_it->x, boundpoint_it->y, boundpoint_it->z;
    }


    if (flag)
        *other_cloud += *ellipse;

    return boundpoint;
}

void path_generater::compute_coverage(Eigen::Vector3d node, double radius){

    std::vector<float> Distance;
    std::vector<int> pointIdx;
    pcl::PointXYZRGB point;
    // float radius = 1;
    point.x = node[0];
    point.y = node[1];
    point.z = node[2];
    kdtree.radiusSearch(point, radius, pointIdx, Distance);
    for(auto& i : pointIdx)
        coverage_flag[i] = 1;

}

// key = 0 up, key = 1 down
int path_generater::compute_boundary(Spline path, Spline* boundary)
{

    double miny = path.miny(), maxy = path.bigy();
    std::map<double, std::vector<double>> boundary_node;
    Eigen::Vector3d point;
    Eigen::Vector3f boundpoint;
    // printf("maxy:%f, miny:%f\n", maxy, miny);
    ////genearte boundary point of path node////
    double dy = miny + 2;
    while(dy < maxy - 2){
        point = path.point(dy);
        dy += toolRadius/4; 
        boundpoint = Area2Cloud(point, 1, 0);
        if(isnan(boundpoint[0])){
            printf("Area Estimatin is NAN\n");
            continue;
        }

        std::vector<double> contain = {boundpoint[0], boundpoint[2]};
        boundary_node[boundpoint[1]] = contain;
    }
    //last point in case in sufficient node
    dy = maxy - 1;
   
    boundpoint = Area2Cloud(point, 1, 0);
    std::vector<double> contain = {boundpoint[0], boundpoint[2]};
    boundary_node[boundpoint[1]] = contain;
    // std::cout << boundpoint[1] << std::endl;
    ////genearte boundary point of path node////

    // spline function of path node
    int node_number = boundary_node.size(), index = 0;
    if(node_number <= 2) return 0;

    double point_x[node_number + 2], point_y[node_number + 2], point_z[node_number + 2];
    for (auto &i : boundary_node){
        index++;
        point_x[index] = i.second[0];
        point_y[index] = i.first;
        point_z[index] = i.second[1]; 
    }
    //get longer boundary, add new node at start&end
    point_x[0] = point_x[1];
    point_y[0] = point_y[1] - 20;
    point_z[0] = point_z[1];
    point_x[index + 1] = point_x[index];
    point_y[index + 1] = point_y[index] + 20;
    point_z[index + 1] = point_z[index];


    Spline bound(node_number + 2, point_y, point_x, point_z);
    *boundary = bound;
    return 1;
}

Eigen::Vector3d path_generater::bisection(Eigen::Vector3d PathNode, Spline* boundary, int itr){

    if(itr > 5) return PathNode;
    Eigen::Vector3f Area_bound; 
    Eigen::Vector3d boundary_point, norm;
    
    //use area_point to fix path node
    //area point is origin_path node's down bounded point
    Area_bound = Area2Cloud(PathNode, 0, 1);
 
    /// check if overflow
    if (Area_bound[1] < boundary->miny() || Area_bound[1] > boundary->bigy())
        return PathNode;
    boundary_point = boundary->point(Area_bound[1]);

    // norm << (Area_bound[0] - boundary_point[0]),
    //     (Area_bound[1] - boundary_point[1]),
    //     (Area_bound[2] - boundary_point[2]);
    norm << (Area_bound[0] - boundary_point[0]), 0, 0;
    if(norm.norm() < Adjust_Threshold) return PathNode;
    PathNode << PathNode[0] - norm[0],
            PathNode[1] - norm[1],
            PathNode[2] - norm[2];

    if(isnan(PathNode[0])){
        std::cout << "adjust path node NAN" << std::endl;
        return PathNode;
    }

    return bisection(PathNode, boundary, itr+1);
}

void path_generater::dynamic_adjust_path(Spline *origin_path, Spline pre_path){

    Spline *boundary = new Spline();
    if(compute_boundary(pre_path, boundary) == 0){
        printf("generate boundary fail\n");
        return;  
    } 
    
    std::map<double, std::vector<double>> new_path;
    double miny = origin_path->miny(), maxy = origin_path->bigy(), dy = 0;
    Eigen::Vector3d path_node;
    std::vector<float> Distance;
    pcl::Indices pointIdx;
    pcl::PointXYZRGB point;

    // printf("origin big y is %f\n", maxy);
    // int NumOfNode = 50;  
    int NumOfNode = (maxy - miny) / 5; 
    for (int i = 1; i < (NumOfNode); i++){

        //注意數量有是否正確
        dy = ((maxy - miny) / NumOfNode * i) + miny;
        path_node = origin_path->point(dy);
        path_node = bisection(path_node, boundary, 0);

        point.x = path_node[0];
        point.y = path_node[1];
        point.z = path_node[2];
        kdtree.nearestKSearch(point, 3, pointIdx, Distance);

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
    }
    // printf("new big y is %f\n", point_y[index - 1]);
    origin_path->restart(index, point_y, point_x, point_z);
    drawpath(*boundary,10,255,10);
    delete boundary;
    boundary = NULL;
}

void path_generater::drawpath(Spline path,int r, int g, int b){

    double miny = path.miny(), maxy = path.bigy(), dy = 0;
    Eigen::Vector3f area_point, com_norm;
    Eigen::Vector3d path_node, bound_point;
    std::vector<float> Distance;
    pcl::Indices pointIdx;
    pcl::PointXYZRGB point;
    int num = 200 ;
    for (int i = 0; i < 200; i++){
        dy = ((maxy - miny) / num * i) + miny;
        path_node = path.point(dy);
        point.x = path_node[0];
        point.y = path_node[1];
        point.z = path_node[2];
        // std::cout << point.x << " " << point.y << std::endl;
        kdtree.nearestKSearch(point, 3, pointIdx, Distance);
        cloud->points[pointIdx[0]].r = (uint8_t)r;
        cloud->points[pointIdx[0]].g = (uint8_t)g;
        cloud->points[pointIdx[0]].b = (uint8_t)b;
    }
}

void path_generater::path_track(Eigen::Vector3f plane_point)
{

    // Insert Point Node with a line
    Eigen::Vector3f normal(1, 0, 0);
    std::map<double, std::vector<double>> Node, boundary;
    Node = insert_point(rangedX_index(plane_point[0]), plane_point);

    int node_number = Node.size(), index = 0;
    double point_x[node_number], point_y[node_number], point_z[node_number];
    for (auto &i : Node)
    {
        point_x[index] = i.second[0];
        point_y[index] = i.first;
        point_z[index] = i.second[1];
        index++;
    }
    Spline path(node_number, point_y, point_x, point_z);
    // drawpath(path, 0, 0, 255);
    // double dy, miny = path.miny(), maxy = path.bigy();
    // Eigen::Vector3d point;

    // for (int i = 0; i < 20; i++)
    // {
    //     dy = ((maxy - miny) / 20 * i) + miny;
    //     point = path.point(dy);
    // }
    Path_set.push_back(path);
}

void path_generater::Contact_Path_Generation()
{
    printf("Start Path Planning!\n");
    struct timeval start, end;
    unsigned long time;
    gettimeofday(&start, NULL);
    kdtree.setInputCloud(cloud);
    //coverage rate
    cloud_size = cloud -> size();
    coverage_flag.resize (cloud_size);
    for(auto &i : coverage_flag){
        i = 0;
    }

    // Plane model
    Eigen::Vector3f normal(1, 0, 0);
    Eigen::Vector3f plane_point(0, 0, 0);

    // find the bounding point
    pcl::PointXYZRGB min_pt, max_pt;
    pcl::getMinMax3D(*cloud, min_pt, max_pt);

    int step_size = toolRadius * 2, path_nums = 0;
    float locateX = min_pt.x + toolRadius;
    // float locateX = 0;
    Spline *boundary = new Spline();
    int GenerateBound = 0;
    while (locateX < max_pt.x){
        plane_point[0] = locateX;
        path_track(plane_point);
        GenerateBound = compute_boundary(Path_set[path_nums], boundary);
        if(GenerateBound == 1) drawpath(*boundary,0,255,0);
        if(locateX != min_pt.x + toolRadius)
            dynamic_adjust_path(&Path_set[path_nums], Path_set[path_nums-1]);
        locateX += step_size;
        drawpath(Path_set[path_nums],255,0,0);
        path_nums++;   
        printf("generate path: %d\n", path_nums);    
    }
    //show surface ande boundary of last one
    // GenerateBound = compute_boundary(Path_set[path_nums-1], boundary);
    // drawpath(*boundary,10,255,10);
    //last straight line path
    // plane_point[0] = max_pt.x - toolRadius / 2;
    // path_track(plane_point);
    // GenerateBound = compute_boundary(Path_set[path_nums], boundary);
    // drawpath(Path_set[path_nums],255,0,0);
    // path_nums++; 
    // printf("generate path(compensation): %d\n", path_nums); 
    delete boundary;

    gettimeofday(&end, NULL);
    time = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    double TIMES = time / 1000000;
    printf("Toal Using Time: %lf\n", TIMES);
    printf("Number of paths: %d\n", path_nums);
}

void path_generater::get_coverage(){

    float yes, no ,rate;
    yes = 0;
    no = 0;
    for(auto i : coverage_flag){
        if(i == 1)
            yes++;
        else
            no++;
    }
    rate = yes / ( yes + no );
    printf("yes: %f, no: %f\n", yes, no);
    printf("coverage rate: %f\n", rate);
}