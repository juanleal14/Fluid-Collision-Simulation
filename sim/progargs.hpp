#ifndef FLUID_OPTIMIZATION_PROGARGS_H
#define FLUID_OPTIMIZATION_PROGARGS_H

#endif //FLUID_OPTIMIZATION_PROGARGS_H
#include "particles.cpp"

static constexpr double const_r = 1.695;
static constexpr double global_density = 1000;
static constexpr double stiff_pressure = 3.0;
static constexpr double stiff_collision = 30000;   //stiffness collision
static constexpr double damping = 128.0;
static constexpr double viscosity = 0.4;
static constexpr double part_size = 0.0002; //Particle Size
static constexpr double time_step = 0.001;
static constexpr double gravity = 9.8;
//std::string address = "../new.fld";





//struct Initial_values;
class Initial_Values;
void check_command_errors(int argc,std::vector<std::string> arguments);
Initial_Values read_general_info(std::ifstream &file);
Particle read_particle(std::ifstream & file);
Grid initialize_grid(std::ifstream &file,Initial_Values &initialValues,int &counter);
Grid initial_read(const std::string& file_address,Initial_Values &initialValues);