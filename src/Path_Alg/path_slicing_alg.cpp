#include "Path_Generate_Algorithm.h"

SectPath::SectPath(std::string configName, std::string CloudFileName)
    :cloud_name(CloudFileName)
{
    read_config(configName);
    cloud = pcl::make_shared<PointCloudType>();
    other_cloud = pcl::make_shared<PointCloudType>();
    normal_cloud = pcl::make_shared<pcl::PointCloud<pcl::Normal>>();
    if (pcl::io::loadPCDFile<pcl::PointXYZRGB>(cloud_name, *cloud) == -1)
        PCL_ERROR("Cloudn't read file!\n");
    else{
        Path_set = {};
        for (size_t i = 0; i < cloud->points.size(); ++i)
        {
            cloud->points[i].r = (uint8_t)255;
            cloud->points[i].g = (uint8_t)255;
            cloud->points[i].b = (uint8_t)255;
            if(ifChangeRange){
                cloud->points[i].x *= 1000;
                cloud->points[i].y *= 1000;
                cloud->points[i].z *= 1000;
            } 
        }
        TransAlign = Eigen::Matrix4f::Identity();
    }
    if(ifSmooth) smooth();
    if(ifAlign) trans2center();
    if(ifRemove) remove_outlier();
};

void SectPath::read_config(std::string filename)
{
    std::ifstream cFile(filename);
    if (!cFile.is_open()) {
        std::cerr << "Couldn't open config file for reading.\n";
        return;
    }

    std::string line;
    while (std::getline(cFile, line)) {
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        auto delimiterPos = line.find("=");
        if (line[0] == '#' || line.empty() || delimiterPos == std::string::npos)
            continue;
        auto name = line.substr(0, delimiterPos);
        auto value = line.substr(delimiterPos + 1);
        if (name == "pathFile") pathFile = value;
        else if (name == "Tool_Radius") toolRadius = std::stod(value);
        else if (name == "PathResolution") PathResolution = std::stod(value);
        else if (name == "RPYresolution") RPYres = std::stod(value);
        else if (name == "Endeffectorlength") EElen = std::stod(value);
        else if (name == "Alignment") {
            if(value == "true") ifAlign = true;
            else  ifAlign = false;
        }else if (name == "Smooth") {
            if(value == "true") ifSmooth = true;
            else ifSmooth = false;
        }else if (name == "ChangeRange") {
            if(value == "true") ifChangeRange = true;
            else ifChangeRange = false;
        }else if (name == "RemoveOutlier") {
            if(value == "true") ifRemove = true;
            else ifRemove = false;
        }     
    }   
}

void SectPath::show()
{
    *other_cloud += *cloud;
    pcl::visualization::PCLVisualizer viewer("VOXELIZED SAMPLES CLOUD");
    viewer.addPointCloud<pcl::PointXYZRGB>(other_cloud, "sample_cloud");

    viewer.setBackgroundColor(0, 0, 0);                                                                        
    viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample_cloud");
    viewer.addCoordinateSystem(100.0);

    viewer.spin();
}

void SectPath::trans2center()
{
    Eigen::Matrix3f covariance_matrix;
    Eigen::Vector4f xyz_centroid4;
    Eigen::Vector3f xyz_centroid3;

    pcl::compute3DCentroid(*cloud, xyz_centroid4);    
    xyz_centroid3 = xyz_centroid4.head<3>();
    pcl::computeCovarianceMatrix(*cloud, xyz_centroid4, covariance_matrix);
   
    Eigen::EigenSolver<Eigen::Matrix3f> eigensolver(covariance_matrix);
    Eigen::VectorXf eigen_values = eigensolver.eigenvalues().real();
    Eigen::Matrix3f eigen_vectors = eigensolver.eigenvectors().real();

    TransAlign.block<3, 1>(0, 3) = -eigen_vectors.transpose()*xyz_centroid3;
    TransAlign.block<3, 3>(0, 0) = eigen_vectors.transpose();
    pcl::transformPointCloud(*cloud, *cloud, TransAlign);
}

void SectPath::remove_outlier()
{
    pcl::StatisticalOutlierRemoval<pcl::PointXYZRGB> sor;
    sor.setInputCloud (cloud);
    sor.setMeanK (50);
    sor.setStddevMulThresh (1.0);
    sor.filter (*cloud);
}


void SectPath::smooth(){

    PointCloudType::Ptr cloud1 = pcl::make_shared<PointCloudType>();
    pcl::search::KdTree<PointType>::Ptr tree = pcl::make_shared<pcl::search::KdTree<PointType>>();
    pcl::PointCloud<pcl::PointNormal> mls_points;
    pcl::MovingLeastSquares<PointType, PointType> mls;
    // mls.setComputeNormals (true);
    // Set parameters
    mls.setInputCloud (cloud);
    mls.setPolynomialOrder (3);
    mls.setSearchMethod (tree);
    mls.setSearchRadius (15);
    // Reconstruct
    mls.setCacheMLSResults(true);
    mls.getCacheMLSResults();
    mls.process(*cloud1);
    pcl::copyPointCloud(*cloud1, *cloud);
    std::string smooth_nam = "smooth_";
    smooth_nam += cloud_name;
    
    if(ifChangeRange){
        for (size_t i = 0; i < cloud1->size(); ++i) {
            cloud1->points[i].x /= 1000;
            cloud1->points[i].y /= 1000;
            cloud1->points[i].z /= 1000;
        }
    }
    pcl::io::savePCDFileASCII(smooth_nam, *cloud1);
}

void SectPath::estimate_normal()
{
    pcl::NormalEstimation<PointType, pcl::Normal> normal_estimation;
    normal_estimation.setInputCloud(cloud);
    pcl::search::KdTree<PointType>::Ptr tree = pcl::make_shared<pcl::search::KdTree<PointType>>();
    normal_estimation.setSearchMethod(tree);
    normal_estimation.setRadiusSearch(2.5);
    // normal_estimation.setKSearch (10);
    normal_estimation.compute(*normal_cloud);
}

std::vector<int> SectPath::rangedX_index(int position)
{
    std::vector<int> indices;
    pcl::PassThrough<PointType> ptfilter(true);
    ptfilter.setInputCloud(cloud);
    ptfilter.setFilterFieldName("x");
    ptfilter.setFilterLimits(-2 + position, 2 + position);
    ptfilter.filter(indices);

    return indices;
}

MAP SectPath::insert_point(std::vector<int> indices, Eigen::Vector3f PlanePoint)
{
    std::vector<int> El, Er;
    std::vector<int> left_pair, right_pair;
    std::vector<int> idr(1), idl(1), idsave(1);
    std::vector<float> dis(1);
    
    // distinguish right or left points of plane
    Eigen::Vector3f point;
    Eigen::Vector3f normal(1, 0, 0);
    for (int &i : indices)
    {
        point[0] = cloud->points[i].x;
        point[1] = cloud->points[i].y;
        point[2] = cloud->points[i].z;
        float distance2plane = (point - PlanePoint).dot(normal);
        if (distance2plane > 0) El.push_back(i);
        else if (distance2plane < 0) Er.push_back(i);
    }
   
    PointCloudType::Ptr cloudEl = pcl::make_shared<PointCloudType>();
    PointCloudType::Ptr cloudEr = pcl::make_shared<PointCloudType>();
    pcl::KdTreeFLANN<PointType> treeEl;
    pcl::KdTreeFLANN<PointType> treeEr;
    pcl::copyPointCloud(*cloud, El, *cloudEl);
    pcl::copyPointCloud(*cloud, Er, *cloudEr);
    treeEl.setInputCloud(cloudEl);
    treeEr.setInputCloud(cloudEr);

    //left set which closet to right, right set which closet to left
    for (int i = 0; i < int(El.size()); i++)
    {
        //left set point
        auto pl = cloudEl->points[i];
        //use left point to search nearest point in right point set
        treeEr.nearestKSearch(pl, 1, idr, dis);
        //right point closeted to left point
        auto pr = cloudEr->points[idr[0]];
        //use right point to search nearest point in left point set
        kdtree.nearestKSearch(pr, 1, idsave, dis);
        right_pair.push_back(idsave[0]);
        
        treeEl.nearestKSearch(pr, 1, idl, dis);
        pl = cloudEl->points[idl[0]];
        kdtree.nearestKSearch(pl, 1, idsave, dis);
        left_pair.push_back(idsave[0]);
    }

    // insert new point to cloud
    PointCloudType insert_cloud;
    insert_cloud.width = left_pair.size();
    insert_cloud.height = 1;
    insert_cloud.points.resize(insert_cloud.width * insert_cloud.height);
    // use for spline
    MAP Node;

    for (int i = 0; i < left_pair.size(); i++)
    {
        int index_right = right_pair[i], index_left = left_pair[i];
        float t = (PlanePoint[0] - cloud->points[index_right].x) / (cloud->points[index_left].x - cloud->points[index_right].x);
        insert_cloud.points[i].x = PlanePoint[0];
        insert_cloud.points[i].y = cloud->points[index_right].y + t * (cloud->points[index_left].y - cloud->points[index_right].y);
        insert_cloud.points[i].z = cloud->points[index_right].z + t * (cloud->points[index_left].z - cloud->points[index_right].z);
        insert_cloud.points[i].r = (uint8_t)255;
        insert_cloud.points[i].g = (uint8_t)0;
        insert_cloud.points[i].b = (uint8_t)0;

        std::vector<double> temp = {insert_cloud.points[i].x, insert_cloud.points[i].z};
        Node[insert_cloud.points[i].y] = temp;
    }
    *other_cloud = (*other_cloud) + insert_cloud;
    
    return Node;
}


Spline SectPath::OnePath(Eigen::Vector3f plane_point)
{
    // Insert Point Node with a line
    MAP Node, boundary;
    int node_number, index = 0;
    double *point_x, *point_y, *point_z;
    
    Node = insert_point(rangedX_index(plane_point[0]), plane_point);
    node_number = Node.size();
    // std::cout << Node.size() << std::endl;

    point_x = new double[node_number];
    point_y = new double[node_number];
    point_z = new double[node_number];
    for (auto &i : Node){
        point_x[index] = i.second[0];
        point_y[index] = i.first;
        point_z[index] = i.second[1];
        index++;
    }
    Spline path(node_number, point_y, point_x, point_z);
    // std::cout << path.miny() << " " << path.bigy() << std::endl;
    delete[] point_x;
    delete[] point_y;
    delete[] point_z;

    return path;
}

void SectPath::drawpath(Spline path,int r, int g, int b)
{
    double miny = path.miny(), maxy = path.bigy();
    std::vector<float> Distance;
    pcl::Indices pointIdx;
    PointType PCLp;

    double dy = miny;
    while(dy < maxy) {
        auto pathPoint = path.point(dy);
        PCLp.x = pathPoint[0];
        PCLp.y = pathPoint[1];
        PCLp.z = pathPoint[2];
        kdtree.nearestKSearch(PCLp, 1, pointIdx, Distance);
        cloud->points[pointIdx[0]].r = (uint8_t)r;
        cloud->points[pointIdx[0]].g = (uint8_t)g;
        cloud->points[pointIdx[0]].b = (uint8_t)b;
        dy += 1; 
    }
}

void SectPath::GenPath()
{
    printf("Start Path Planning!\n");
    auto startT = std::chrono::high_resolution_clock::now();

    kdtree.setInputCloud(cloud);
    //coverage rate
    // coverage_flag.resize (cloud -> size());
    // std::fill(coverage_flag.begin(), coverage_flag.end(), 0);
    // Plane model
    Eigen::Vector3f plane_point(0, 0, 0);
    // find the bounding point
    PointType min_pt, max_pt;
    pcl::getMinMax3D(*cloud, min_pt, max_pt);
    std::cout << min_pt.y << " " << max_pt.y << std::endl;
    ///////////////////////////////////////////////////////
    // >>>>>>>>>>>>>>>>  generate path  <<<<<<<<<<<<<<<< //
    std::vector<Spline> path_set2front, path_set2back;
    int step_size = toolRadius * 2, path_nums = 0;
    // generate front path //
    float loc = (min_pt.x + max_pt.x)/2 - step_size;
    while (loc > min_pt.x){
        plane_point[0] = loc;
        Spline loc_path = OnePath(plane_point);
        path_set2front.insert(path_set2front.begin(), loc_path);
        drawpath(loc_path, 255, 0, 0);
        loc -= step_size;
        path_nums++;   
        printf("generate path: %d\n", path_nums);   ;  
    }
    // generate front path //
    loc = (min_pt.x + max_pt.x)/2;
    while (loc < max_pt.x){
        plane_point[0] = loc;
        Spline loc_path = OnePath(plane_point);
        path_set2back.push_back(loc_path);
        drawpath(loc_path, 255, 0, 0);
        loc += step_size;
        path_nums++;   
        printf("generate path: %d\n", path_nums);    
    }
    Path_set = std::move(path_set2front);
    std::move(path_set2back.begin(), path_set2back.end(), std::back_inserter(Path_set));
    // >>>>>>>>>>>>>>>>  generate path  <<<<<<<<<<<<<<<< //
    ///////////////////////////////////////////////////////


    auto endT = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endT - startT).count();
    printf("Toal Using Time: %ld \n", duration);
    printf("Number of paths: %d\n", path_nums);
    printf("Number of Point Cloud: %ld\n", cloud -> size());
}
