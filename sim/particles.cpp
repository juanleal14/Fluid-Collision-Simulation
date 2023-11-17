#include "particles.hpp"
struct Particle{
    double px;
    double py;
    double pz;
    double hvx;
    double hvy;
    double hvz;
    double vx;
    double vy;
    double vz;
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

