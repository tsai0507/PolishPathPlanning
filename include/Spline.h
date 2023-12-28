#ifndef SPLINE
#define SPLINE

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

class Spline{
public:
    Spline(){};
    Spline(int number, const double* point_y, const double* point_x, const double* point_z): 
        node_number(number), big_y(point_y[number-1]), small_y(point_y[0])
    {   
        acc = gsl_interp_accel_alloc ();
        Method = gsl_interp_steffen;
        splineYX = gsl_spline_alloc (Method, node_number);
        splineYZ = gsl_spline_alloc (Method, node_number);
        gsl_spline_init (splineYX, point_y, point_x, node_number);
        gsl_spline_init (splineYZ, point_y, point_z, node_number);

    };
    //must use double type (Vector3d)
    Eigen::Vector3d point(double y){
        Eigen::Vector3d value(gsl_spline_eval(splineYX, y, acc), y, gsl_spline_eval(splineYZ, y, acc));
        return value;
    }

    double miny() { return small_y;}
    double bigy() { return big_y;}

    void restart(int number, const double* point_y, const double* point_x, const double* point_z){       
        gsl_spline_free(splineYX);
        gsl_spline_free(splineYZ);
        gsl_interp_accel_free (acc);
        acc = gsl_interp_accel_alloc();
        node_number = number;
        splineYX = gsl_spline_alloc (Method, node_number);
        splineYZ = gsl_spline_alloc (Method, node_number);
        gsl_spline_init (splineYX, point_y, point_x, node_number);
        gsl_spline_init (splineYZ, point_y, point_z, node_number);
        big_y = point_y[number-1];
        small_y = point_y[0];
    }

private:
    int node_number;
    double big_y, small_y;
    gsl_interp_accel *acc;
    const gsl_interp_type *Method;
    gsl_spline *splineYX;
    gsl_spline *splineYZ;
};

class PPP{
public:
    PPP(){};
    PPP(int number, double* point_y): 
        node_number(number)
    {   
        point = point_y[0];
    };


    double miny() { return node_number;}
    double bigy() { return node_number;}

private:
    int node_number;
    double point;
};

#endif