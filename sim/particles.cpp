#include "particles.hpp"

 template <typename T>
 class Vect3
 {
public:
    // Constructor to initialize constant features
    Vect3(T xset, T yset, T zset) : coords({xset,yset,zset})  {}
    // Getter methods to access individual features
    [[nodiscard]] T x() const { return coords[0]; }
    [[nodiscard]] T y() const { return coords[1]; }
    [[nodiscard]] T z() const { return coords[2]; }
    void set_x(T x)  { coords[0] = x; }
    void set_y(T y)  { coords[1] = y; }
    void set_z(T z)  { coords[2] = z; }
    void set(T x, T y, T z) { coords[0] = x; coords[1] = y; coords[2] = z; }
    T& operator[](int i){return coords[i];}
    Vect3 operator- (Vect3<T>& other) {
        Vect3 result(coords[0]-other[0],coords[1]-other[1],coords[2]-other[2]);
        return result;
    }

    Vect3 operator* (const double &c) const{
        Vect3 const result(coords[0]*c,coords[1]*c,coords[2]*c);
        return result;
    }

    Vect3 operator= (Vect3<T> &other){ ///Operator =
        std::copy(coords.begin(),coords.end(),&other[0]);
        return *this;
    }
    bool operator== (Vect3<T> &other){ ///Operator =
        bool result = (coords[0]==other[0]) && (coords[1]==other[1]) && (coords[2]==other[2]);
        return result;
    }
    double dist_sqrd(Vect3& other){
        const Vect3<T> diff = *this - other;
        return pow(diff.x(),2) + pow(diff.y(),2) + pow(diff.z(),2);

    }


private:
    std::array<T,3> coords; ///Checkear si se aloca junto a otros parametros de la clase o separado
};


class Particle{
public:
        int id;
    Vect3<double> pos;
    Vect3<double> hv;
    Vect3<double> v;
    double density;
    Vect3<double> acceleration;

    bool operator== (Particle &other){ ///Operator =
        bool cond = (id == other.id) && (pos == other.pos) && (hv == other.hv) && (v == other.v) && (density == other.density) && (acceleration == other.acceleration);
        return cond;
    }

    Particle() : id(0), pos(0.0, 0.0, 0.0), hv(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0), density(0.0), acceleration(0.0, 0.0, 0.0) {};
    Particle(int id_val, Vect3<double> pos, Vect3<double> hvs, Vect3<double> vel, double d, Vect3<double> a): id(id_val), pos(pos), hv(hvs), v(vel), density(d), acceleration(a){}
    double distance_to(Particle& p2){return pos.dist_sqrd(p2.pos);

    /*
    Vect3 pos(){return position_values;}
    Vect3 hv(){return hv_values;}
    Vect3 vel(){return v_values;}
    int id(){return id_number;}
    Vect3 acc()
     {return acceleration;}
     */

};



