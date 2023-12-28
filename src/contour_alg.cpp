#include "contour_alg.h"

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
    // cloud = cloud1;
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
        insert_cloud.points[i].g = (uint8_t)255;
        insert_cloud.points[i].b = (uint8_t)255;

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
        printf("generate path: %d\n", path_nums);    
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
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endT - startT).count()*0.001;
    printf("Toal Using Time: %lf (ms)\n", duration);
    printf("Number of paths: %d\n", path_nums);
    printf("Number of Point Cloud: %ld\n", cloud -> size());
}

void SectPath::HandEyeTransform(Eigen::VectorXf &point)
{
    // >>>>>>>>>>>>>>> HandEye matrix <<<<<<<<<<<<<<< //
    Eigen::Vector3f position(HANDEYEx, HANDEYEy, HANDEYEz);
    Eigen::Matrix4f HEmatrix = Eigen::Matrix4f::Identity();
    Eigen::Matrix3f HErot;
    //// Euler rotation Z-Y-X ////
    HErot = Eigen::AngleAxisf(HANDEYErz, Eigen::Vector3f::UnitZ())* 
            Eigen::AngleAxisf(HANDEYEry, Eigen::Vector3f::UnitY())*
            Eigen::AngleAxisf(HANDEYErx, Eigen::Vector3f::UnitX());
    HEmatrix.block<3, 3>(0, 0) = HErot;
    HEmatrix.block<3, 1>(0, 3) = position;

    // >>>>>>>>>>>>>>> coordinate matrix <<<<<<<<<<<<<<< //
    Eigen::Matrix4f coordinate = Eigen::Matrix4f::Identity();
    Eigen::Matrix3f posture;
    //// Euler rotation Z-Y-X ////
    posture = Eigen::AngleAxisf(point[5], Eigen::Vector3f::UnitZ())* 
                Eigen::AngleAxisf(point[4], Eigen::Vector3f::UnitY())*
                Eigen::AngleAxisf(point[3], Eigen::Vector3f::UnitX());
    coordinate.block<3, 3>(0, 0) = posture;
    coordinate.block<3, 1>(0, 3) = point.head<3>();

    // >>>>>>>>>>>>>>> Coordinate Transformation <<<<<<<<<<<<<<< //
    coordinate = HEmatrix*coordinate;
    Eigen::Vector3f translation = coordinate.block<3, 1>(0, 3);
    Eigen::Matrix3f rotation = coordinate.block<3, 3>(0, 0);
    Eigen::Vector3f euler_angles = rotation.eulerAngles(2, 1, 0);

    point << translation.x(), translation.y(), translation.z(),
            euler_angles(2), euler_angles(1), euler_angles(0);
    // std::cout << coordinate << std::endl;
    // std::cout << "---------------------------" << std::endl;
}

void SectPath::reduceRPY(std::vector<int> Index)
{
    if(RPYres <= 2) return;

    int preId = 0, lastId = 0, keyId = 0, res = RPYres;
    for(int id = 0; id < Index.size(); id++){  
        do{
            double dr[3];
            lastId = preId + res;
            for(size_t D = 3; D < 6; D++){
                if(WayPointsList[lastId][D]*WayPointsList[preId][D] >= 0){
                    dr[D-3] = (WayPointsList[lastId][D] - WayPointsList[preId][D])/res;
                }else{
                    double no1, no2;
                    if(WayPointsList[lastId][D] < 0){
                        no2 = WayPointsList[preId][D];
                        no1 = 2*M_PI + WayPointsList[lastId][D];
                    }else{
                        no2 = 2*M_PI + WayPointsList[preId][D];
                        no1 = WayPointsList[lastId][D];
                    }
                    dr[D-3] = std::abs(WayPointsList[lastId][D] - WayPointsList[preId][D]) < std::abs(no1-no2)
                                ?(WayPointsList[lastId][D] - WayPointsList[preId][D]) :(no1-no2);
                    dr[D-3] /= res;
                }    
            }
            for(int wi = 1; wi < res; wi++){
                WayPointsList[preId+wi][3] = dr[0] + WayPointsList[preId+wi-1][3];
                WayPointsList[preId+wi][4] = dr[1] + WayPointsList[preId+wi-1][4];
                WayPointsList[preId+wi][5] = dr[2] + WayPointsList[preId+wi-1][5];
            }
            preId = lastId;
        }while((preId+res) <= Index[id]);

        if(preId != Index[id]){
            for(int i = preId+1; i <= Index[id]; i++){
                WayPointsList[i][3] = WayPointsList[preId][3];
                WayPointsList[i][4] = WayPointsList[preId][4];
                WayPointsList[i][5] = WayPointsList[preId][5];
            }
        }
        preId = Index[id] + 1; 
    }
    //limit RPY in -180~180
    for(auto &point: WayPointsList){    
        point[3] = point[3]>M_PI ?point[3]-2*M_PI :point[3];
        point[4] = point[4]>M_PI ?point[4]-2*M_PI :point[4];
        point[5] = point[5]>M_PI ?point[5]-2*M_PI :point[5];
    }
}


void SectPath::TransFlangeposition()
{
    //// TransEE2Flange //// 
    Eigen::Vector3f EE2Flange(0, 0, -EElen);
    Eigen::Matrix4f EE2Fmatrix = Eigen::Matrix4f::Identity();
    EE2Fmatrix.block<3, 1>(0, 3) = EE2Flange.head<3>();
    
    for(auto &point: WayPointsList){
        Eigen::Matrix4f coordinate = Eigen::Matrix4f::Identity();
        Eigen::Matrix3f posture;
        //// Euler rotation Z-Y-X ////
        posture = Eigen::AngleAxisf(point[5], Eigen::Vector3f::UnitZ())* 
                    Eigen::AngleAxisf(point[4], Eigen::Vector3f::UnitY())*
                    Eigen::AngleAxisf(point[3], Eigen::Vector3f::UnitX());
        coordinate.block<3, 3>(0, 0) = posture;
        coordinate.block<3, 1>(0, 3) = point.head<3>();
        coordinate = coordinate*EE2Fmatrix;

        Eigen::Vector3f translation = coordinate.block<3, 1>(0, 3);
        point[0] = translation.x();
        point[1] = translation.y();
        point[2] = translation.z();
    }
}

void SectPath::postion_smooth()
{
    /*** Path positon(xyz) smooth ***/
    std::vector<Eigen::VectorXf> new_path(WayPointsList);
    double tolerance = 0.00001, change = tolerance, weight_data = 0.65, weight_smooth = 1-weight_data;
    double x_i, y_i, y_prev, y_next, y_i_saved;
    int dim = 3, path_len = WayPointsList.size();

  	while(change >= tolerance)
    {
  		change = 0;
  		for(int i = 1; i < path_len - 1; i++){
  			for(int j = 0; j < dim; j++){
                x_i = WayPointsList[i][j];
                y_i = new_path[i][j];
                y_prev = new_path[i - 1][j];
                y_next = new_path[i + 1][j];

                y_i_saved = y_i;
                y_i += (weight_data * (x_i - y_i) + weight_smooth * (y_next + y_prev - 2 * y_i));

                new_path[i][j] = y_i;
  				change += abs(y_i - y_i_saved);
  			}
  		}
  	}
    WayPointsList = std::move(new_path);
}


void SectPath::getPath()
{
    Eigen::Matrix4f invTransAlign = TransAlign.inverse();
    std::vector<std::vector<Eigen::Vector4f>> WayPointsXYZList;
    int flag = 1;
    Path_set.erase(Path_set.begin());
    Path_set.pop_back();
    std::cout << "Path number is " << Path_set.size() << std::endl;
    // std::for_each(Path_set.begin(), Path_set.end(), [](Spline path) 
    //     { std::cout << path.miny() << std::endl; });

    for(auto path: Path_set){
        std::vector<Eigen::Vector4f> oneWayPathXYZ;
        double dy = path.miny() + 5;
        while(dy < path.bigy() - 5){
            auto pointXYZ = path.point(dy);
            Eigen::Vector4f wayPointXYZ(pointXYZ[0], pointXYZ[1], pointXYZ[2], 1);
            wayPointXYZ = invTransAlign * wayPointXYZ;
            oneWayPathXYZ.push_back(wayPointXYZ);
            dy += PathResolution; 
        }
        if(flag == -1) std::reverse(oneWayPathXYZ.begin(), oneWayPathXYZ.end());
        WayPointsXYZList.push_back(oneWayPathXYZ);
        flag *= -1;
    }
    // *cloud += *other_cloud;
    pcl::transformPointCloud(*cloud, *cloud, invTransAlign);
    pcl::transformPointCloud(*other_cloud, *other_cloud, invTransAlign);
    kdtree.setInputCloud(cloud);
    this->estimate_normal();

    // compute RPY and output final path
    std::vector<int> TailIndex; //use for reduceRPY
    for(auto onePathXYZ: WayPointsXYZList){
        for(auto PointXYZ: onePathXYZ){
            std::vector<int> id;
            std::vector<float> dis;
            float raw, pitch, yall;
            pcl::PointXYZRGB SPoint;
            SPoint.x = PointXYZ[0];
            SPoint.y = PointXYZ[1];
            SPoint.z = PointXYZ[2];
            Eigen::VectorXf WayPoint(6);

            kdtree.nearestKSearch(SPoint, 1, id, dis);
            pcl::Normal N(normal_cloud->points[id[0]]);

            Eigen::Vector3f Approach(float(-N.normal_x), float(-N.normal_y), float(-N.normal_z));
            Eigen::Vector3f Orientation = Approach.cross(Eigen::Vector3f::UnitX());;
            Eigen::Vector3f Normal = Orientation.cross(Approach);
            Eigen::Matrix3f rotationMatrix = Eigen::Matrix3f::Identity();
            rotationMatrix.col(0) = Normal;
            rotationMatrix.col(1) = Orientation;
            rotationMatrix.col(2) = Approach;
            Eigen::Vector3f euler_angles = rotationMatrix.eulerAngles(2, 1, 0);
            raw = euler_angles(2);
            pitch = euler_angles(1);
            yall = euler_angles(0);

            if(ifChangeRange)
                WayPoint << PointXYZ[0]/1000, PointXYZ[1]/1000, PointXYZ[2]/1000, raw, pitch, yall;
            else
                WayPoint << PointXYZ[0], PointXYZ[1], PointXYZ[2], raw, pitch, yall;
            HandEyeTransform(WayPoint);
            WayPointsList.push_back(WayPoint);
        } TailIndex.push_back(WayPointsList.size()-1);
    }
    postion_smooth();
    reduceRPY(TailIndex);
    TransFlangeposition();

    std::cout << "!!!!! GOT PATH !!!!!" << std::endl;
    std::ofstream outputFile(pathFile);
    if (outputFile.is_open()){
        for (const auto& point : WayPointsList){
            for (int i = 0; i < 6; i++)
                outputFile << point[i] << " ";
            outputFile << std::endl;
        }
        outputFile.close();
        std::cout << "File saved: " << pathFile << std::endl;
    }
    else
        std::cerr << "Unable to open file: " << pathFile << std::endl;
}




