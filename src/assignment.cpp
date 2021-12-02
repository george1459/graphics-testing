#include "object.h"
#include "scene.h"

#include <iostream>

#include "image.h"

using namespace Eigen;
using namespace std;

const int MAX_ITERS = 10000;
const int XRES = 500;
const int YRES = 500;

/**
 * IOTest Code
 */

bool Superquadric::IOTest(const Vector3d &point) {
    Matrix4d O = Matrix4d::Identity();
    for (auto &i: transforms) {
        O = i.get()->GetMatrix() * O;
    }
    Vector4d new_point;
    new_point << (point)[0], (point)[1], (point)[2], 1.0;
    new_point = O.inverse() * new_point;

    double x = new_point[0];
    double y = new_point[1];
    double z = new_point[2];

    double S = -1.0 + pow(z*z, 1.0/exp1) + pow(pow(x*x, 1.0/exp0) + pow(y*y, 1.0/exp0), exp0/exp1);
    // cout << S << endl;

    if (S < 0) {
        return true;
    }

    /**
     * PART 1
     * TODO: Implement the IO Test function for a superquadric. Remember to
     *       apply any transformations to the superquadric before testing
     *       against the IO function.
     */
    return false;
}

bool Assembly::IOTest(const Vector3d &point) {
    Matrix4d O = Matrix4d::Identity();
    for (auto &i: transforms) {
        O = i.get()->GetMatrix() * O;
    }
    Vector4d new_point;
    new_point << (point)[0], (point)[1], (point)[2], 1.0;
    new_point = O.inverse() * new_point;

    double x = new_point[0];
    double y = new_point[1];
    double z = new_point[2];

    Vector3d new_point_input;
    new_point_input << x, y, z;

    bool res = false;
    for (auto &i: children) {
        res = res || i.get()->IOTest(new_point_input);
    }
    /**
     * PART 1
     * TODO: Implement the IO Test function for an assembly (recursively call
     *       IOTest on the children). Make sure to apply any transformations
     *       to the assembly before calling IOTest on the children.
     */
    return res;
}

double S_func(double x, double y, double z, double e, double n) {
    return -1.0 + pow(z*z, 1.0/n) + pow(pow(x*x, 1.0/e) + pow(y*y, 1.0/e), e/n);
}

double S_func_t(double ax, double ay, double az, double bx, double by, double bz,
double t, double e, double n) {
    return S_func(ax*t+bx, ay*t+by, az*t+bz, e, n);
}

double S_deriv_func(double ax, double ay, double az, double bx, double by, double bz,
double t, double e, double n) {
    double x_comp = ax * t + bx;
    double y_comp = ay * t + by;
    double z_comp = az * t + bz;
    // double res = 2*az*bz*pow(pow(z_comp,2.0),-1.0+1.0/n) + 2*az*az*t*pow(pow(z_comp,2.0),-1.0+1.0/n)
    // + (ax*pow(x_comp,-1.0+1.0/e)+ay*pow(y_comp,-1.0+1.0/e))*pow(pow(x_comp,1.0/e)+pow(y_comp,1.0/e),-1.0+e/n);
    // return res/n;
    double common = pow(pow(x_comp*x_comp, 1.0/e) + pow(y_comp*y_comp, 1.0/e), e/n - 1.0);
    Vector3d later;
    later << 2.0 * x_comp * pow(x_comp*x_comp, 1.0/e - 1.0) * common, 2.0 * y_comp * pow(y_comp*y_comp, 1.0/e - 1) * common, 2.0 * z_comp * pow(z_comp*z_comp, 1.0/n - 1.0);
    Vector3d a_vec;
    a_vec << ax, ay, az;
    return a_vec.dot(later)/n;
}

/**
 * Closest Intersection Code
 */

pair<double, Intersection> Superquadric::ClosestIntersection(const Ray &ray) {

    Matrix4d O = Matrix4d::Identity();
    for (auto &i: transforms) {
        O = i.get()->GetMatrix() * O;
    }
    Vector4d new_point;
    Vector3d a_vec, b_vec;
    new_point << (ray.origin)[0], (ray.origin)[1], (ray.origin)[2], 1.0;
    new_point = O.inverse() * new_point;
    b_vec << new_point[0], new_point[1], new_point[2];

    new_point << (ray.direction)[0], (ray.direction)[1], (ray.direction)[2], 1.0;
    new_point = O.inverse() * new_point;
    a_vec << new_point[0], new_point[1], new_point[2];

    // Derive the a, b, c parameters as in lecture notes
    double a, b, c;
    a = a_vec.dot(a_vec);
    b = 2 * a_vec.dot(b_vec);
    c = b_vec.dot(b_vec) - 3;
    int sign_b;
    sign_b = (b >= 0) ? 1 : -1;

    if (b*b - 4*a*c < 0) {
        pair<double, Intersection> closest = make_pair(INFINITY, Intersection());
        return closest;
    }

    // get t_1 and t_2
    double t_1, t_2, t;
    t_1 = (-b - sign_b * sqrt(b*b - 4*a*c))/(2*a);
    t_2 = (2*c)/(-b - sign_b * sqrt(b*b - 4*a*c));

    // using the smaller value if it is positive
    double initial_t = min(t_1, t_2);

    // TODO: maybe check this later
    if (initial_t < 0) {
        pair<double, Intersection> closest = make_pair(INFINITY, Intersection());
        return closest;
    }
    
    // initial_t = (t_1 < t_2 && t_1 > 0)? t_1: t_2;
    t = initial_t;

    double n = exp1;
    double e = exp0;

    double ax, ay, az, bx, by, bz;
    ax = a_vec[0];
    ay = a_vec[1];
    az = a_vec[2];
    bx = b_vec[0];
    by = b_vec[1];
    bz = b_vec[2];

    while (abs(S_func_t(ax,ay,az,bx,by,bz,t,e,n)) >= 1.0/12.0) {
        double deriv = S_deriv_func(ax,ay,az,bx,by,bz,t,e,n);
        double value = S_func_t(ax,ay,az,bx,by,bz,t,e,n);
        if (deriv > 0) {
            break;
        }
        if (abs(deriv) < 1E-6) {
            // "if the value of the function is no where close to 0"
            if (abs(value) > 1E-6) {
                pair<double, Intersection> closest = make_pair(INFINITY, Intersection());
                return closest;
            }
            break;
        }
        t = t - value/deriv;
    }

    /**
     * PART 1
     * TODO: Implement a ray-superquadric intersection using Newton's method.
     *       Make sure to apply any transformations to the superquadric before
     *       performing Newton's method.
     */

    // Now do calculations for the intersection part
    Vector3d intersection_point;
    intersection_point << ray.direction[0] * t + ray.origin[0], ray.direction[1] * t + ray.origin[1], ray.direction[2] * t + ray.origin[2];

    Vector3d normal, transformed_intersection;
    Vector4d temp;
    temp << intersection_point[0], intersection_point[1], intersection_point[2], 1;
    temp = O.inverse() * temp;
    transformed_intersection << temp[0], temp[1], temp[2];

    normal = GetNormal(transformed_intersection);
    normal = O.inverse().transpose().block<3, 3>(0, 0) * normal;

    Ray new_ray;
    new_ray.direction = normal;
    new_ray.origin = transformed_intersection;

    pair<double, Intersection> closest = make_pair(t, Intersection(new_ray, this));
    return closest;
}

pair<double, Intersection> Assembly::ClosestIntersection(const Ray &ray) {
    /**
     * PART 1
     * TODO: Implement a ray-assembly intersection by recursively finding
     *       intersection with the assembly's children. Make sure to apply any
     *       transformations to the assembly before calling ClosestIntersection
     *       on the children.
     */
    pair<double, Intersection> closest = make_pair(INFINITY, Intersection());
    for (auto &i: children) {
        pair<double, Intersection> res = i.get()->ClosestIntersection(ray);
        if (res.first <= closest.first) {
            closest = res;
        }
    }

    return closest;
}

/**
 * Raytracing Code
 */

void Scene::Raytrace() {
    Image img = Image(XRES, YRES);

    for (int i = 0; i < XRES; i++) {
        for (int j = 0; j < YRES; j++) {
            /**
             * PART 2
             * TODO: Implement raytracing using the code from the first part
             *       of the assignment. Set the correct color for each pixel
             *       here.
             */
            img.SetPixel(i, j, Vector3f::Ones());
        }
    }

    // Outputs the image.
    if (!img.SaveImage("rt.png")) {
        cerr << "Error: couldn't save PNG image" << std::endl;
    } else {
        cout << "Done!\n";
    }
}
