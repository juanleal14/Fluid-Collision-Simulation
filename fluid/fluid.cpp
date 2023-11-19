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
    //simulate(1,grid,initialValues);
    //write_to_file(arguments[2],grid,initialValues);

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

    load_trace("./trz/small/denstransf-base-1.trz",grid_trz,initialValues);
    //for (int block = 0; block< grid_trz.blocks.size(); block++){for (auto part: grid_trz[block]){if (part.id==2635){std::cout<<"Bloque"<<block<<'\n';}}}
    accelerations_transfer(grid_trz,initialValues);
    check_trace("./trz/small/acctransf-base-1.trz",grid_trz);
    /*
    load_trace("./trz/small/motion-base-1.trz",grid_trz,initialValues);
    Vect3<int> belongings (0,0,0);
    for (int i = 0; i<grid_trz.blocks.size(); i++){
      belongings = belongs_to_boundary_block(i,grid_trz.size);
      for (int j = 0; j<grid_trz[i].size(); j++) {
        new_boundary_interaction(belongings,grid_trz[i][j]);
     }
    }
    check_trace("./trz/small/boundint-base-1.trz",grid_trz);
    */

    //check_trace("./trz/small/acctransf-base-1.trz",grid_trz);
    //densities_increase(grid,initialValues);
    //check_trace("./trz/small/densinc-base-1.trz",grid);

}
