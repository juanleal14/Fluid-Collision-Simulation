#include "libraries.h"
#include "grid.hpp"


void compare_accelerations(Vect3 &a1, Vect3 &a2, long id);
void compare_particle(Particle &p1, Particle &p2,long id);
void find_elem(int e, std::vector<int> &v);
void check_trace(std::string trz, Grid &grid, std::vector<Particle> &particles, std::vector<double> &densities, std::vector<Acceleration> &accelerations);
std::vector <Particle> receive_trace(std::string trz, std::vector<double> &densities, std::vector<Acceleration> &accelerations);
void write_to_file(const std::string& output_file_address,std::vector<Particle> particles, Initial_Values &initialValues);
void particles_motion(std::vector<Particle> &particles, Grid &grid, std::vector <Acceleration> &accelerations);
void acceleration_transfer(std::vector<Particle> &particles, Grid &grid, std::vector<double> &densities, std::vector<Acceleration> &accelerations);
void accelerations_computation(std::vector<Particle> &particles, Grid &grid,std::vector <double> &densities, std::vector <Acceleration> &accelerations);
void densities_transform(std::vector<double> &densities);
void densities_increase(std::vector<Particle> &particles, Grid &grid, std::vector<double> &densities);

void simulate(int nsteps, std::vector<Particle> &particles, Grid &grid);
