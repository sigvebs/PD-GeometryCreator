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
    int q = 20;
    int threshold = 2000;
    int imageResolution = 6000;

    // Boundaries
    double X_0 = 0;
    double X_1 = 1;
    double Y_0 = 0;
    double Y_1 = 1;

    // alpha beta
    double alpha_1 = 0.5, alpha_2 = 0.5;
    double beta_1 = 0.5, beta_2 = 0.5;

    bool periodic_x = false;
    bool periodic_y = false;

    string basePath = "/media/Media4/Scratch/MeshGenerator/tmp";
    string imgPath = "";

    bool testingSave = false;
    int testSaveFreq = 100;

    int redistributionFrequency = 100;
    int nRedistributedPoints = 0;
};
//------------------------------------------------------------------------------
class MeshGenerator
{
public:
    MeshGenerator();
    MeshGenerator(Parameters parameters);
    void initializeFromImage();
    arma::mat createMesh();

    void createDomainGrid();
    void mapParticlesToGrid();
    void save_image_and_xyz(string base, int nr = -1);
    void setDomainSize(double spacing);
    double calculateRadialDistribution(int nr = -1);
    void writeConfiguration();
protected:
    Parameters param;

    // From image
    CImg<double> image;

    int h;
    int w;
    arma::mat img_data;

    int n;
    int q;
    int threshold;

    double alpha_1;
    double alpha_2;
    double beta_1;
    double beta_2;
    double X_0;
    double X_1;
    double Y_0;
    double Y_1;

    arma::mat x;
    arma::vec js;
    std::vector<vector<int>> gridNeighbours;
    std::vector<vector<int>> particlesInGridPoint;

    unsigned seed;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution_x;
    std::uniform_real_distribution<double> distribution_y;


    double dx;
    double dy;
    double DX;
    double DY;

    bool periodic_x;
    bool periodic_y;

    double optimalGridSpacing;

    // Domain variables
    int dim = 2;
    double nx, ny;
    double gridSpacing_x, gridSpacing_y;

    string basePath;
    int imageResolution;
    int testSaveFreq;

    int findGridId(const arma::vec2 & r_i);
    void checkBoundaries();

    int openmp_threads;
};
//------------------------------------------------------------------------------
// Functions
//------------------------------------------------------------------------------
//void save_xyz(arma::mat &x, std::string base, int i);
//void save_image_and_xyz(Parameters parameters, arma::mat &x, std::string base, int i);
//void save_image_and_xyz(int resolution, arma::mat &image, arma::mat &x, std::string base, int i);
//arma::mat createMeshFromImage(Parameters parameters);
//------------------------------------------------------------------------------
}
#endif // MESHGENERATOR_H
