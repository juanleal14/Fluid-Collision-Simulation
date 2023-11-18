#ifndef FLUID_PARTICLES_MOTION_H
#define FLUID_PARTICLES_MOTION_H
#endif
#include "grid.cpp"


//void compare_accelerations(Vect3<double> &a1, Vect3<double> &a2, long id);
//void compare_particle(Particle &p1, Particle &p2,long id);
//void find_elem(long id, Block &v);
//void check_trace(std::string trz, Grid &grid);
//std::vector <Particle> receive_trace(std::string trz, Grid &grid);
void write_to_file(const std::string &output_file_address,Grid grid, Initial_Values &initialValues);
void particles_motion(Grid &grid);
void acceleration_transfer(Grid &grid);
void accelerations_computation(Grid &grid);
void densities_transform(Grid &grid);
void densities_increase(Grid &grid);
void simulate(int nsteps, Grid &grid);

