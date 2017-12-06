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
    try {
        cfg.readFile( cfgFileName.c_str() );
    }
    catch (const libconfig::FileIOException &fioex) {
        std::cerr << "I/O error while reading the configuration file." << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (const libconfig::ParseException &pex) {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        exit(EXIT_FAILURE);
    }
    libconfig::Setting &root = cfg.getRoot();

    //--------------------------------------------------------------------------
    struct mg::Parameters param;

    // Reading parameters

    if(root.exists("dataType"))
         root.lookupValue("dataType", param.dataType);
    else {
        std::cerr << "'Save path'dataType' not supplied" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(root.exists("nParticles"))
        param.nParticles = root["nParticles"];
    else
        std::cerr << "Number of particles not supplied"<< std::endl;

    if(root.exists("steps"))
        param.threshold = root["steps"];
    else
        std::cerr << "Threshold not supplied" << std::endl;

    if(root.exists("samples"))
        param.q = root["samples"];
    else
        std::cerr << "MultiplicationFactor not supplied" << std::endl;

    if(root.exists("savePath"))
         root.lookupValue("savePath", param.savePath);
    else {
        std::cerr << "Save path not supplied" << std::endl;
        exit(EXIT_FAILURE);
    }

    if(root.exists("dataPath"))
        root.lookupValue("dataPath", param.dataPath);
    else {
        std::cerr << "'dataPath' not supplied" << std::endl;
        exit(EXIT_FAILURE);
    }


    if(root.exists("periodic_x"))
        param.periodic_x = (int) root["periodic_x"];
    if(root.exists("periodic_y"))
        param.periodic_y = (int) root["periodic_y"];
    if(root.exists("periodic_z"))
        param.periodic_z = (int) root["periodic_z"];
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
    if(root.exists("openmp_threads"))
        param.openmp_threads = root["openmp_threads"];
    if(root.exists("imageResolution"))
        param.imageResolution = root["imageResolution"];
    if(root.exists("imageFormat"))
        root.lookupValue("imageFormat", param.imageFormat);

    if(param.dataType == "binvox") { }

    if(param.dataType == "images") {
        if(root.exists("leading_zeros"))
            param.leading_zeros = root["leading_zeros"];

        if(root.exists("z_scale"))
            param.z_scale = root["z_scale"];

        if(root.exists("baseName"))
            root.lookupValue("baseName", param._baseName);
        else {
            std::cerr << "'baseName' not supplied" << std::endl;
            exit(EXIT_FAILURE);
        }

        if(root.exists("voxelRange")) {
            const libconfig::Setting & vRange = root["voxelRange"];
            param.imgRange[0] = vRange[0];
            param.imgRange[1] = vRange[1];
        }else{
            std::cerr << "'voxelRange' not supplied" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(param.dataType == "image") {
        if(root.exists("z_depth"))
            root.lookupValue("z_depth", param.z_depth);
        else{
            param.z_depth = 1;
            std::cout << "z_depth set to 1" << std::endl;
        }
    }

    param.q = param.nParticles*param.q;
    //--------------------------------------------------------------------------

    cout << param.savePath << endl;
    boost::filesystem::path dir(param.savePath);
    if(boost::filesystem::create_directories(dir)) {
        std::cout << "Directory created: " << param.savePath << "\n";
    }

    arma::wall_clock timer;

    timer.tic();

    mg::MeshGenerator mg(param);
    std::cout << "MeshGenerator object generated" << std::endl;
    arma::mat x = mg.createMesh();
    std::cout << "Geometry created" << std::endl;
    std::cout << "Calculating Radial Distribution" << std::endl;
    mg.calculateRadialDistribution();
    std::cout << "Saving data" << std::endl;
    mg.save_image_and_xyz(param.savePath + "/mesh");
    std::cout << "Writing configuration" << std::endl;
    mg.writeConfiguration();

    double n_secs = timer.toc();

    std::cout << "Geometry computed in " << n_secs << " seconds" << std::endl;
    std::cout << "Complete" << std::endl;
    return EXIT_SUCCESS;
}
//------------------------------------------------------------------------------
