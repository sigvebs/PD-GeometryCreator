#ifndef MESHGENERATOR_H
#define MESHGENERATOR_H

#include <stdio.h>
#include <armadillo>
#include <vector>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <random>
#include <chrono>
#include <omp.h>

#include "mg_functions.h"

#include <CImg.h>
using namespace cimg_library;
using namespace std;

//------------------------------------------------------------------------------
// NAMESPACE MG
//------------------------------------------------------------------------------
namespace mg
{
//------------------------------------------------------------------------------
struct Parameters
{
    int nParticles = 128;
    int q = 1024;
    int threshold = 2000;
    int imageResolution = 1024;

    // Boundaries
    double X_0 = 0;
    double X_1 = 1;
    double Y_0 = 0;
    double Y_1 = 1;

    // alpha beta
    double alpha_1 = 0.5, alpha_2 = 0.5;
    double beta_1 = 0.5, beta_2 = 0.5;

    // Radius center hole
    double r_centerHole = 0.2;
    bool periodic_x = false;
    bool periodic_y = false;

    string basePath = "/media/Media4/Scratch/MeshGenerator/tmp";
    string imgPath = "";
};
//------------------------------------------------------------------------------
class MeshGenerator
{
public:
    MeshGenerator();
protected:
    Parameters param;
};
//------------------------------------------------------------------------------
// Functions
//------------------------------------------------------------------------------
void save_xyz(arma::mat &x, std::string base, int i);
void save_image_and_xyz(Parameters parameters, arma::mat &x, std::string base, int i);
void save_image_and_xyz(Parameters parameters, arma::mat &image, arma::mat &x, std::string base, int i);
arma::mat createMeshFromImage(Parameters parameters);
//------------------------------------------------------------------------------
}
#endif // MESHGENERATOR_H
