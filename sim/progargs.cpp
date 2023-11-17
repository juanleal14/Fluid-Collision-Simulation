#include "progargs.hpp"

class Initial_Values {
public:
    // Constructors
    Initial_Values() : ppm(0.0f), np(0), m(0.0), h(0.0) {}

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



void check_command_errors(int argc,char** argv) {
    if (argc != 4) {
        std::cerr << "Error: Invalid number of arguments: " << argc << ".";
        exit(-1);
    }
    //Is integer
    try {
        // Convert the input parameter to an integer
        int value = std::stoi(argv[1]);
        // Use the integer value as needed
        std::cout << "The entered integer is: " << value << std::endl;
        if (value < 0) {
            std::cerr << "Error: The entered integer is negative." << std::endl;
            exit(-2);
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: The entered parameter is not an integer." << std::endl;
        exit(-1);
    }
}


Initial_Values read_general_info(std::ifstream &file){
    // Read ppm and np
    //cap 10-11
    Initial_Values initialValues;
    double read_value;
    file.read(reinterpret_cast<char*>(&read_value), sizeof(float));//NOLINT
    initialValues.setPpm(read_value);
    file.read(reinterpret_cast<char*>(&read_value), sizeof(int));//NOLINT
    initialValues.setNp(read_value);
    initialValues.setM(global_density/pow(initialValues.getPpm(),3));
    initialValues.setH(const_r/initialValues.getPpm());
    std::cout << "ppm: " << initialValues.getPpm() << ", np: " << initialValues.getNp() << std::endl;
    return initialValues;
}
int read_particle_info(std::ifstream &file,std::vector<Particle> &particles){
    int counter = 0;
    //file.read(reinterpret_cast<char*)
    float part_coord_x = 0;
    float part_coord_y = 0;
    float part_coord_z = 0;
    float part_vel_x = 0;
    float part_vel_y = 0;
    float part_vel_z = 0;
    float hvx = 0;
    float hvy = 0;
    float hvz = 0;
    while(file.peek()!=EOF){
        counter ++;
        Particle particle;
        file.read(reinterpret_cast<char*>(&part_coord_x), sizeof(float));//NOLINT
        file.read(reinterpret_cast<char*>(&part_coord_y), sizeof(float));//NOLINT
        file.read(reinterpret_cast<char*>(&part_coord_z), sizeof(float));//NOLINT
        file.read(reinterpret_cast<char*>(&hvx), sizeof(float));//NOLINT
        file.read(reinterpret_cast<char*>(&hvy), sizeof(float));//NOLINT
        file.read(reinterpret_cast<char*>(&hvz), sizeof(float));//NOLINT
        file.read(reinterpret_cast<char*>(&part_vel_x), sizeof(float));//NOLINT
        file.read(reinterpret_cast<char*>(&part_vel_y), sizeof(float));//NOLINT
        file.read(reinterpret_cast<char*>(&part_vel_z), sizeof(float));//NOLINT
        particle.pos.set_x(static_cast<double>(part_coord_x));
        particle.pos.set_y(static_cast<double>(part_coord_y));
        particle.pos.set_z(static_cast<double>(part_coord_z));
        particle.hv.set_x(static_cast<double>(hvx));
        particle.hv.set_y(static_cast<double>(hvy));
        particle.hv.set_z(static_cast<double>(hvz));
        particle.v.set_x(static_cast<double>(part_vel_x));
        particle.v.set_y(static_cast<double>(part_vel_y));
        particle.v.set_z(static_cast<double>(part_vel_z));
        particles.push_back(particle);
    }
    return counter;}
std::vector<Particle> initial_read(const std::string& file_address,Initial_Values &initialValues){
    //Read file
    std::ifstream file(file_address, std::ios::binary);
    if (!file.is_open()) { //Check error opening
        std::cerr<<"Error: Cannot open " << file_address <<" for reading";
        exit (-3);
    }
    std::vector<Particle> particles; //Iniatilize a vector of particles to be stored
    initialValues = read_general_info(file);//call to a function to read parameters
    std::cout<<"Lee bien ppm np";
    int const counter = read_particle_info(file,particles);
    if(counter == 0){
        std::cout<< "Error : Invalid number of particles: " << counter <<".";
    }
    else if(counter != initialValues.getNp()){
        std::cout<<"Error : Number of particles mismatch. Header " << initialValues.getNp() << " Found " << counter <<".";
    }
    file.close();
    return particles;
}


