#include "particles.hpp"

 template <typename T>
 class Vect3
 {
public:
    // Constructor to initialize constant features
    Vect3(T xset, T yset, T zset) : xcoord(xset),ycoord(yset), zcoord(zset)  {}
    // Getter methods to access individual features
    [[nodiscard]] T x() const { return xcoord; }
    [[nodiscard]] T y() const { return ycoord; }
    [[nodiscard]] T z() const { return zcoord; }
    void set_x(T x)  { xcoord = x; }
    void set_y(T y)  { ycoord = y; }
    void set_z(T z)  { zcoord = z; }

    Vect3 operator- (const Vect3 &other) const{
        Vect3 const result(xcoord-other.x(),ycoord-other.y(),zcoord-other.z());
        return result;
    }

    Vect3 operator* (const double &c) const{
        Vect3 const result(xcoord*c,ycoord*c,zcoord*c);
        return result;
    }

    Vect3& operator= (Vect3 const &other){
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
    Vect3<double> pos;
    Vect3<double> hv;
    Vect3<double> v;
    double density;
    Vect3<double> acceleration;

public:
    Particle() : id(0), pos(0.0, 0.0, 0.0), hv(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0), density(0.0), acceleration(0.0, 0.0, 0.0) {};
    Particle(int id_val, Vect3<double> pos, Vect3<double> hvs, Vect3<double> vel, double d, Vect3<double> a): id(id_val), pos(pos), hv(hvs), v(vel), density(d), acceleration(a){}
    double distance_to(Particle& p2){return pos.dist_sqrd(p2.pos);}
    /*
    Vect3 pos(){return position_values;}
    Vect3 hv(){return hv_values;}
    Vect3 vel(){return v_values;}
    int id(){return id_number;}
    Vect3 acc()
     {return acceleration;}
     */
};



