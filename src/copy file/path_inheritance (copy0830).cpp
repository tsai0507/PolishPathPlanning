#include "path_connect.h"

path_generater::path_generater(std::string configName, std::string CloudFileName)
{
    cloud_name = CloudFileName;
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

void path_generater::read_config(std::string filename)
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
        else if (name == "depth") depth = std::stod(value);
        else if (name == "Adjust_Threshold") Adjust_Threshold = std::stod(value);
        else if (name == "toolthickness") toolthickness = std::stod(value);
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

Eigen::Matrix4f path_generater::compute_transform(PointType point, float *principle_curvature)
{

    Eigen::Vector3f curvatureVector, normalVector;
    std::vector<float> Distance;
    pcl::Indices pointIdx;
    pcl::PrincipalCurvaturesEstimation
        <PointType, pcl::Normal, pcl::PrincipalCurvatures> principal_curvatures_estimation;

    // kdtree.radiusSearch(point, toolRadius*0.5, pointIdx, Distance);
    kdtree.nearestKSearch(point, 10, pointIdx, Distance);
    principal_curvatures_estimation.computePointPrincipalCurvatures(
        *normal_cloud, pointIdx[0], pointIdx, curvatureVector[0], curvatureVector[1],
        curvatureVector[2], principle_curvature[0], principle_curvature[1]);

    normalVector << normal_cloud->points[pointIdx[0]].normal_x,
        normal_cloud->points[pointIdx[0]].normal_y,
        normal_cloud->points[pointIdx[0]].normal_z;
    Eigen::Vector3f crossVector = normalVector.cross(curvatureVector);

    // noa Rotation Matrix
    // bigger curvature dir, smallerr curvature dir, normal, point
    Eigen::Matrix4f Y, X(Eigen::Matrix<float, 4, 4>::Identity());
    Y << crossVector(0), curvatureVector(0), normalVector(0), point.x,
        crossVector(1), curvatureVector(1), normalVector(1), point.y,
        crossVector(2), curvatureVector(2), normalVector(2), point.z,
        0, 0, 0, 1;
    Eigen::Matrix4f pt2Base = X.transpose() * Y;

    return pt2Base;
}

// flag draw the data, key up/down
Eigen::Vector3f path_generater::Area2Cloud(Eigen::Vector3d point, bool flag, Dir key)
{
    Eigen::Matrix4f transform;
    float principle_curvature[2];
    double longAxis, shortAxis;
    PointType SearchPoint;

    SearchPoint.x = point[0];
    SearchPoint.y = point[1];
    SearchPoint.z = point[2];
    transform = compute_transform(SearchPoint, principle_curvature);
    
    
    if((principle_curvature[0] >= 0) && (principle_curvature[1] >= 0)){
        longAxis = std::sqrt(
            std::pow(1 / principle_curvature[1], 2) - std::pow(std::abs(1 / principle_curvature[1]) - depth, 2));
        if (longAxis > toolRadius) longAxis = toolRadius;

        shortAxis = std::sqrt(
            std::pow(1 / principle_curvature[0], 2) - std::pow(std::abs(1 / principle_curvature[0]) - depth, 2));
        if (shortAxis > toolRadius) shortAxis = toolRadius;
        
    }
    else{
        longAxis = std::abs(
            1 / principle_curvature[1]) - std::sqrt(std::pow(1 / principle_curvature[1], 2) - std::pow(toolRadius, 2));
        if (longAxis > toolthickness) longAxis = toolthickness;
        
        shortAxis = std::abs(
            1 / principle_curvature[0]) - std::sqrt(std::pow(1 / principle_curvature[0], 2) - std::pow(toolRadius, 2));
        if (shortAxis > toolthickness) shortAxis = toolthickness;
    }
    
    // draw a Ellipse in origin ,and trasform to path node
    auto ellipse = std::make_shared<PointCloudType>();
    std::uint8_t r(255), g(0), b(0);
    for (float angle(0.0); angle <= 360.0; angle += 0.5)
    {
        PointType basic_point;
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
    if (key == right) {
        auto boundpoint_it = std::max_element(ellipse->begin(), ellipse->end(), [](auto &p1, auto &p2)
                                              { return p1.x < p2.x; });
        boundpoint << boundpoint_it->x, boundpoint_it->y, boundpoint_it->z;

        //use for compute coverage
        // auto for_min = std::min_element(ellipse->begin(), ellipse->end(), [](auto &p1, auto &p2)
        //                                       { return p1.x < p2.x; });
        // double comput_lan =   (for_min-> x - boundpoint_it->x)/2; 
        // compute_coverage(point, comput_lan);                           
    }
    else {
        auto boundpoint_it = std::min_element(ellipse->begin(), ellipse->end(), 
            [](auto &p1, auto &p2){ return p1.x < p2.x; });
        boundpoint << boundpoint_it->x, boundpoint_it->y, boundpoint_it->z;
    }
    if(flag) *other_cloud += *ellipse;

    return boundpoint;
}

// key = 0 up, key = 1 down
int path_generater::compute_boundary(Spline path, std::shared_ptr<Spline>& boundary, Dir key)
{
    
    double miny = path.miny(), maxy = path.bigy();
    MAP boundary_node;
    Eigen::Vector3d point;
    Eigen::Vector3f boundpoint;

    ////genearte boundary point of path node////
    double dy = miny + 2;
    while(dy < maxy - 2){
        point = path.point(dy);
        dy += toolRadius/4; 
        boundpoint = Area2Cloud(point, 1, key);
        if(isnan(boundpoint[0])){
            printf("Area Estimatin is NAN\n");
            continue;
        }
        std::vector<double> contain = {boundpoint[0], boundpoint[2]};
        boundary_node[boundpoint[1]] = contain;
    }

    //last point in case in sufficient node
    dy = maxy - 1;
    boundpoint = Area2Cloud(point, 1, key);
    std::vector<double> contain = {boundpoint[0], boundpoint[2]};
    boundary_node[boundpoint[1]] = contain;

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
    
    
    boundary.reset(new Spline(node_number + 2, point_y, point_x, point_z));
    
    return 1;
}

Eigen::Vector3d path_generater::bisection(
    Eigen::Vector3d PathNode, std::shared_ptr<Spline> boundary, int itr, Dir key)
{
    if(itr > 5) return PathNode;
    Eigen::Vector3f Area_bound; 
    Eigen::Vector3d boundary_point, norm;
    
    //use area_point to fix path node
    //area point is origin_path node's down bounded point
    Area_bound = Area2Cloud(PathNode, 0, key == left? right: left);
 
    /// check if overflow
    if (Area_bound[1] < boundary->miny() || Area_bound[1] > boundary->bigy())
        return PathNode;
    boundary_point = boundary->point(Area_bound[1]);

    norm << (Area_bound[0] - boundary_point[0]), 0, 0;
    if(norm.norm() < Adjust_Threshold) return PathNode;
    PathNode << PathNode[0] - norm[0],
            PathNode[1] - norm[1],
            PathNode[2] - norm[2];

    if(isnan(PathNode[0])){
        std::cout << "adjust path node NAN" << std::endl;
        return PathNode;
    }

    return bisection(PathNode, boundary, itr+1, key);
}

void path_generater::dynamic_adjust_path(
        Spline *origin_path, std::shared_ptr<Spline> pre_boundary, Dir key)
{   
    std::map<double, std::vector<double>> new_path;
    double miny = origin_path->miny(), maxy = origin_path->bigy(), dy = 0;
    Eigen::Vector3d path_node;
    std::vector<float> Distance;
    pcl::Indices pointIdx;
    pcl::PointXYZRGB point;

    // int NumOfNode = 50;  
    int NumOfNode = (maxy - miny) / 5; 
    
    
    for (int i = 0; i <= NumOfNode; i++)
    {
        //注意數量有是否正確
        dy = ((maxy - miny) / NumOfNode * i) + miny;
        path_node = origin_path->point(dy);
        path_node = bisection(path_node, pre_boundary, 0, key);

        point.x = path_node[0];
        point.y = path_node[1];
        point.z = path_node[2];
        kdtree.nearestKSearch(point, 3, pointIdx, Distance);

        std::vector<double> contain = {cloud->points[pointIdx[0]].x, cloud->points[pointIdx[0]].z};        
        new_path[cloud->points[pointIdx[0]].y] = contain;
    }
 
    int index = 0;
    double point_x[new_path.size()], point_y[new_path.size()], point_z[new_path.size()];
    for(auto& i: new_path){
        point_x[index] = i.second[0];
        point_y[index] = i.first;
        point_z[index] = i.second[1];
        index++;
    }
    origin_path->restart(index, point_y, point_x, point_z);
}

void path_generater::thread_worker(thread_wrap_data data, std::vector<Spline> &_path_set)
{
    int step_size = toolRadius*2, loc = (data.max_pt + data.min_pt)/2;
    auto boundary = std::make_shared<Spline>();
    Spline pre_path = data.CenPath;
    loc += data.key == left? -step_size: step_size;

    while (data.max_pt > loc && loc > data.min_pt)
    {
        mtx.lock();
        /// critical section which will coredump ///
        Spline loc_path = OnePath(Eigen::Vector3f(loc, 0, 0));
        if(compute_boundary(pre_path, boundary, data.key)) 
            drawpath(*boundary, 0,255,0);
        /// critical section which will coredump ///
        mtx.unlock();
        dynamic_adjust_path(&loc_path, boundary, data.key);

        _path_set.push_back(loc_path);
        pre_path = loc_path;
        drawpath(loc_path, 0, 0, 255);
        if(data.key == left) loc -= step_size;
        else loc += step_size;  
    }
}


void path_generater::GenPath()
{
    printf("Start Path Planning!\n");
    auto startT = std::chrono::high_resolution_clock::now();

    PointType min_pt, max_pt;
    estimate_normal();
    kdtree.setInputCloud(cloud);
    pcl::getMinMax3D(*cloud, min_pt, max_pt); // find the bounding point

    ///////////////////////////////////////////////////////
    // >>>>>>>>>>>>>>>>  generate path  <<<<<<<<<<<<<<<< //
    std::vector<Spline> path_set2front, path_set2back;
    Spline Center_path;

    // generate center path //
    Center_path = OnePath(Eigen::Vector3f((min_pt.x+max_pt.x)/2, 0, 0));
    drawpath(Center_path, 0, 0, 255);

    // thread data
    thread_wrap_data data1(left, min_pt.x, max_pt.x, Center_path);
    thread_wrap_data data2(right, min_pt.x, max_pt.x, Center_path);
    // thread work
    std::thread thread1(std::bind(
        &path_generater::thread_worker, this, data1, std::ref(path_set2front)));
    std::thread thread2(std::bind(
        &path_generater::thread_worker, this, data2, std::ref(path_set2back)));
    // wait thread 
    thread1.join();
    thread2.join();

    std::move(path_set2front.rbegin(), path_set2front.rend(), std::back_inserter(Path_set));
    Path_set.push_back(std::move(Center_path));
    std::move(path_set2back.begin(), path_set2back.end(), std::back_inserter(Path_set));
    // >>>>>>>>>>>>>>>>  generate path  <<<<<<<<<<<<<<<< //
    ///////////////////////////////////////////////////////

    auto endT = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endT - startT).count()*0.001;
    printf("Toal Using Time: %lf (ms)\n", duration);
    // printf("Number of paths: %d\n", path_nums);
    printf("Number of Point Cloud: %ld\n", cloud -> size());
}