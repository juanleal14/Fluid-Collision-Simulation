#ifndef FLUID_OPTIMIZATION_PROGARGS_H
#define FLUID_OPTIMIZATION_PROGARGS_H

#endif //FLUID_OPTIMIZATION_PROGARGS_H
#include "particles.cpp"



//std::string address = "../new.fld";





//struct Initial_values;
class Initial_Values;
void check_command_errors(int argc,std::vector<std::string> arguments);
Initial_Values read_general_info(std::ifstream &file);
Particle read_particle(std::ifstream & file);
