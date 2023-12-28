#include "Path_Generate_Algorithm.h"

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
    // Path_set.pop_back();
    std::cout << "Path number is " << Path_set.size() << std::endl;
    // std::for_each(Path_set.begin(), Path_set.end(), [](Spline path) 
    //     { std::cout << path.miny() << std::endl; });

    for(auto path: Path_set){
        std::vector<Eigen::Vector4f> oneWayPathXYZ;
        double dy = path.miny() + 10;
        while(dy < path.bigy() - 10){
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