#include "sim/particles_motion.cpp"
//using namespace std;



int main(int argc, char** argv) {
    check_command_errors(argc,argv);
    Initial_Values initialValues;
    //const std::string text = argv[2];
    std::span const args_view{argv, static_cast<std::size_t>(argc)};
    std::vector<std::string> const arguments{args_view.begin() + 1, args_view.end()};
    check_command_errors(argc,arguments);
    Grid grid = initial_read(arguments[2],initialValues);
    std::cout<<"\nNum particles: "<<grid.size()<<'\n';
    Grid grid = grid_initialization(initialValues);
    simulate(1,grid);
    write_to_file(arguments[3],grid,initialValues);

// TRAZE  BOUNDARY ITERATION
    std::vector<double> densities;
    std::vector<Acceleration> accelerations;

    load_trace("./trz/small/acctransf-base-1.trz",grid,myparticles,densities,accelerations,initialValues);
    particle_collision(myparticles,grid,accelerations);
    check_trace("./trz/small/partcol-base-1.trz",grid,myparticles,densities,accelerations);
    //check_trace("./trz/small/boundint-base-1.trz",grid,myparticles,densities,accelerations);
}
