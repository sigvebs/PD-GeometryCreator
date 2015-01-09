#include <iostream>
#include <armadillo>
#include <boost/filesystem.hpp>
#include <libconfig.h++>

#include "../src/meshgenerator.h"
using namespace std;

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    libconfig::Config cfg;

    string cfgFileName;
    if (argc == 2){
        cfgFileName = argv[1];
    }else{
        std::cerr << "Supply a valid path for the configuration file" << std::endl;
        return EXIT_FAILURE;
    }

    // Reading configuration file
    try
    {
        cfg.readFile( cfgFileName.c_str() );
    }
    catch (const libconfig::FileIOException &fioex)
    {
        std::cerr << "I/O error while reading the configuration file." << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (const libconfig::ParseException &pex)
    {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        exit(EXIT_FAILURE);
    }
    libconfig::Setting &root = cfg.getRoot();

    //--------------------------------------------------------------------------
    struct mg::Parameters param;

    // Reading parameters

    if(root.exists("nParticles"))
        param.nParticles = root["nParticles"];
    else
        std::cerr << "Number of particles not supplied"<< std::endl;

    if(root.exists("threshold"))
        param.threshold = root["threshold"];
    else
        std::cerr << "Threshold not supplied" << std::endl;

    if(root.exists("multiplicationFactor"))
        param.q = root["multiplicationFactor"];
    else
        std::cerr << "MultiplicationFactor not supplied" << std::endl;

    if(root.exists("savePath"))
         root.lookupValue("savePath", param.basePath);
    else{
        std::cerr << "Save path not supplied" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(root.exists("imgPath"))
        root.lookupValue("imgPath", param.imgPath);
    else{
        std::cerr << "Image path not supplied" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(root.exists("periodic_x"))
        param.periodic_x = root["periodic_x"];
    if(root.exists("periodic_y"))
        param.periodic_y = root["periodic_y"];
    if(root.exists("alpha_1"))
        param.alpha_1 = root["alpha_1"];
    if(root.exists("alpha_2"))
        param.alpha_2 = root["alpha_2"];
    if(root.exists("beta_1"))
        param.beta_1 = root["beta_1"];
    if(root.exists("beta_2"))
        param.beta_2 = root["beta_2"];
    if(root.exists("debug"))
        param.testingSave = root["debug"];
    if(root.exists("testSaveFreq"))
        param.testSaveFreq = root["testSaveFreq"];
    if(root.exists("redistributionFrequency"))
        param.redistributionFrequency = root["redistributionFrequency"];
    if(root.exists("nRedistributedPoints"))
        param.nRedistributedPoints = root["nRedistributedPoints"];

    param.q = param.nParticles*param.q;
    //--------------------------------------------------------------------------

    boost::filesystem::path dir(param.basePath);
    if(boost::filesystem::create_directories(dir)) {
        std::cout << "Directory created: " << param.basePath << "\n";
    }

    arma::wall_clock timer;

    timer.tic();

    mg::MeshGenerator mg(param);
    arma::mat x = mg.createMesh();
    std::cout << "Geometry created" << std::endl;
    mg.save_image_and_xyz(param.basePath + "/mesh");
    mg.calculateRadialDistribution();
    mg.writeConfiguration();

    double n_secs = timer.toc();

    std::cout << "Geometry computed in " << n_secs << " seconds" << std::endl;
    std::cout << "Complete" << std::endl;
    return EXIT_SUCCESS;
}
//------------------------------------------------------------------------------
