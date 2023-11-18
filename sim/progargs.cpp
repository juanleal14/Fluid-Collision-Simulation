#include "progargs.hpp"
#include "grid.cpp"


class Initial_Values {
public:
    // Constructors
    Initial_Values() : ppm(0.0), np(0), m(0.0), h(0.0) {}

    Initial_Values(float ppmValue, int npValue, double mValue, double hValue)
            : ppm(ppmValue), np(npValue), m(mValue), h(hValue) {}

    // Getter methods
    [[nodiscard]] float getPpm() const { return ppm; }
    [[nodiscard]] int getNp() const { return np; }
    [[nodiscard]] double getM() const { return m; }
    [[nodiscard]] double getH() const { return h; }

    // Setter methods
    void setPpm(float ppmValue) { ppm = ppmValue; }
    void setNp(int npValue) { np = npValue; }
    void setM(double mValue) { m = mValue; }
    void setH(double hValue) { h = hValue; }
private:
    float ppm;
    int np;
    double m;
    double h;
};



void check_command_errors(int argc,std::vector<std::string> arguments) {
    if (argc != 4) {
        std::cerr << "Error: Invalid number of arguments: " << argc << ".";
        exit(-1);
    }
    //Is integer
    try {
        // Convert the input parameter to an integer
        int const value = std::stoi(arguments[1]);
        // Use the integer value as needed
        std::cout << "The entered integer is: " << value << "\n";
        if (value < 0) {
            std::cerr << "Error: The entered integer is negative.\n";
            exit(-2);
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: The entered parameter is not an integer.\n";
        exit(-1);
    }
}


Initial_Values read_general_info(std::ifstream &file){
    // Read ppm and np
    //cap 10-11
    Initial_Values initialValues;
    float read_value;
    file.read(reinterpret_cast<char*>(&read_value), sizeof(float));//NOLINT
    initialValues.setPpm(read_value);
    int read_value2;
    file.read(reinterpret_cast<char*>(&read_value2), sizeof(int));//NOLINT
    initialValues.setNp(read_value2);
    initialValues.setM(global_density/pow(initialValues.getPpm(),3));
    initialValues.setH(const_r/initialValues.getPpm());
    std::cout << "ppm: " << initialValues.getPpm() << ", np: " << initialValues.getNp() << std::endl;
    return initialValues;
}

Particle read_particle(std::ifstream & file){ ///Read only one particle from file
    float value = 0;
    Particle particle;
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.pos.set_x(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.pos.set_y(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.pos.set_z(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.hv.set_x(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.hv.set_y(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.hv.set_z(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.v.set_x(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.v.set_y(static_cast<double>(value));
    file.read(reinterpret_cast<char*>(&value), sizeof(float));//NOLINT
    particle.v.set_z(static_cast<double>(value));

    return particle;
}

Grid initialize_grid(std::ifstream &file,Initial_Values &initialValues,int &counter){ /// Crear Grid con el size correcto y añadir todas las partículas del archivo

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
    Grid grid(gridSize);

    while(file.peek()!=EOF){
        counter ++;
        grid.add_particle(read_particle(file));
    }
    std::cout<<"Total particles read from file = "<<counter<<'\n';
    return grid;
}


Grid initial_read(const std::string& file_address,Initial_Values &initialValues){
    //Read file
    std::ifstream file(file_address, std::ios::binary);
    int counter = 0;
    if (!file.is_open()) { //Check error opening
        std::cerr<<"Error: Cannot open " << file_address <<" for reading";
        exit (-3);
    }
    initialValues = read_general_info(file);//call to a function to read parameters
    std::cout<<"Lee bien ppm np";
    Grid grid = initialize_grid(file,initialValues,counter);
    if(counter == 0){
        std::cout<< "Error : Invalid number of particles: " << counter <<".";
    }
    else if(counter != initialValues.getNp()){
        std::cout<<"Error : Number of particles mismatch. Header " << initialValues.getNp() << " Found " << counter <<".";
    }
    file.close();
    return grid;
}


