#include "sim/particles_motion.cpp"
//using namespace std;



int main(int argc, char** argv) {
    check_command_errors(argc,argv);
    Initial_Values initialValues;
    const std::string text = argv[2];
    std::vector<Particle> myparticles = initial_read(argv[2],initialValues);
    std::cout<<"\nNum particles: "<<myparticles.size()<<'\n';
    //h = initialValues.getH();
    //m = initialValues.getM();
    Grid grid = grid_initialization(initialValues,myparticles);
    simulate(1,myparticles,grid);
    //write_to_file(argv[3],myparticles,initialValues);

// TRAZE  BOUNDARY ITERATION
    std::vector<double> densities;
    std::vector<Acceleration> accelerations;

    load_trace("./trz/small/acctransf-base-1.trz",grid,myparticles,densities,accelerations,initialValues);
    particle_collision(myparticles,grid,accelerations);
    check_trace("./trz/small/partcol-base-1.trz",grid,myparticles,densities,accelerations);
    //check_trace("./trz/small/boundint-base-1.trz",grid,myparticles,densities,accelerations);
}
