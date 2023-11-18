#include "sim/particles_motion.cpp"

//using namespace std;



int main(int argc, char** argv) {
    Initial_Values initialValues;
    //const std::string text = argv[2];
    std::span const args_view{argv, static_cast<std::size_t>(argc)};
    std::vector<std::string> const arguments{args_view.begin() + 1, args_view.end()};
    check_command_errors(argc,arguments);
    //std::cout<<"\nNum particles: "<<grid.size.<<'\n';
    Grid grid = initial_read(arguments[1],initialValues);
    simulate(1,grid,initialValues);
    write_to_file(arguments[2],grid,initialValues);

    const double boxx = bmax_coord_x - bmin_coord_x;
    const double boxy = bmax_coord_y - bmin_coord_y;
    const double boxz = bmax_coord_z - bmin_coord_z;
    GridSize gridSize;
    gridSize.setNumX(floor(boxx/initialValues.getH()));
    gridSize.setNumY(floor(boxy/initialValues.getH()));
    gridSize.setNumZ(floor(boxz/initialValues.getH()));
    gridSize.setSizeX(boxx/gridSize.getNumX());
    gridSize.setSizeY(boxy/gridSize.getNumY());
    gridSize.setSizeZ(boxz/gridSize.getNumZ());
    Grid grid_trz(gridSize) ;

    load_trace("./trz/small/acctransf-base-1.trz",grid_trz,initialValues);
    particle_collision(grid_trz);
    check_trace("./trz/small/partcol-base-1.trz",grid_trz);
    //check_trace("./trz/small/boundint-base-1.trz",grid,myparticles,densities,accelerations);
}
