#ifndef FLUID_OPTIMIZATION_PROGARGS_H
#define FLUID_OPTIMIZATION_PROGARGS_H

#endif //FLUID_OPTIMIZATION_PROGARGS_H
#include "particles.cpp"

constexpr double const_r = 1.695;
constexpr double global_density = 1000;
constexpr double stiff_pressure = 3.0;
constexpr double stiff_collision = 30000;   //stiffness collision
constexpr double damping = 128.0;
constexpr double viscosity = 0.4;
constexpr double part_size = 0.0002; //Particle Size
constexpr double time_step = 0.001;
constexpr double gravity = 9.8;
//std::string address = "../new.fld";





//struct Initial_values;
class Initial_Values;
void check_command_errors(int argc,char** argv);
Initial_Values read_general_info(std::ifstream &file);
Particle read_particle(std::ifstream & file);
