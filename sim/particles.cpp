#include "particles.hpp"

class Vect3 {
public:
    // Constructor to initialize constant features
    Vect3(double xset, double yset, double zset) : xcoord(xset), zcoord(zset), ycoord(yset) {}
    // Getter methods to access individual features
    [[nodiscard]] T x() const { return xcoord; }
    [[nodiscard]] T y() const { return ycoord; }
    [[nodiscard]] T z() const { return zcoord; }
    void set_x(T x)  { x = x; }
    void set_y(T y)  { y = y; }
    void set_z(T z)  { z = z; }

    Vect3 operator- (const Vect3 &other) const{
        Vect3 const result(xcoord-other.x(),ycoord-other.y(),zcoord-other.z());
        return result;
    }

    Vect3 operator* (const double &c) const{
        Vect3 const result(xcoord*c,ycoord*c,zcoord*c);
        return result;
    }

    Vect3& operator= (Vect3& other){
        xcoord = other.x();
        ycoord = other.y();
        zcoord = other.z();

        return *this;
    }

    double dist_sqrd(Vect3& other){
        Vect3 diff = *this - other;
        return pow(diff.x(),2) + pow(diff.y(),2) + pow(diff.z(),2);

    }


private:
    double xcoord{};
    double ycoord{};
    double zcoord{};
};


class Particle{
public:
    int id;
    Vect3 pos;
    Vect3 hv;
    Vect3 v;
    double density;
    Vect3 acceleration;

public:
    Particle(int id_val, Vect3 pos, Vect3 hvs, Vect3 vel, double d, Vect3 a): id(id_val), pos(pos), hv(hvs), v(vel), density(d), acceleration(a){}
    /*
    Vect3 pos(){return position_values;}
    Vect3 hv(){return hv_values;}
    Vect3 vel(){return v_values;}
    int id(){return id_number;}
    Vect3 acc(){return acceleration;}
     */
};

struct Acceleration{
    double ax;
    double ay;
    double az;
};
double distance_squared(Particle p1, Particle p2) {
    double dx = p1.px - p2.px;
    double dy = p1.py - p2.py;
    double dz = p1.pz - p2.pz;

    return dx * dx + dy * dy + dz * dz;
}

