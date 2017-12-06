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
#include "grid.h"

#include <CImg.h>
using namespace cimg_library;
using namespace std;

//------------------------------------------------------------------------------
// NAMESPACE MG
//------------------------------------------------------------------------------
namespace mg
{
    typedef arma::Cube<unsigned char> cube_type;
//static const int NORM_N = 2;
#define NORM_N 2
//------------------------------------------------------------------------------
struct Parameters
{
    int nParticles = 128;
    int q = 20;
    int threshold = 2000;
    int imageResolution = 6000;
    int z_depth = 1;

    // Boundaries
    double X_0 = 0;
    double X_1 = 1;
    double Y_0 = 0;
    double Y_1 = 1;
    double Z_0 = 0;
    double Z_1 = 1;

    // alpha beta
    double alpha_1 = 0.5, alpha_2 = 0.5;
    double beta_1 = 0.5, beta_2 = 0.5;

    bool periodic_x = false;
    bool periodic_y = false;
    bool periodic_z = false;

    string dataPath = "/media/Media4/Scratch/MeshGenerator/tmp";

    bool testingSave = false;
    int testSaveFreq = 100;
    string savePath;

    int redistributionFrequency = 100;
    int nRedistributedPoints = 0;

    int openmp_threads = 2;

    string _baseName;
    int imgRange[2];
    double z_scale = 1.0;
    int leading_zeros = 3;

    string dataType = "binvox";
    string imageFormat = "png";

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
    void setDomainSize(double spacing, double solidFraction);
    double calculateRadialDistribution(int nr = -1);
    void writeConfiguration();

    void readBinvox();
    void readImage();
    void readImages();
    void printProgress(double progress);
protected:
    void computeGridPorosity();
    void computePorosity(Grid & grid, const cube_type &img_data);

    Parameters param;

    int m_h;
    int m_w;
    int m_d;
    cube_type m_img_data;
    double m_solidFraction;

    int m_nParticles;
    int m_samples;
    int m_steps;

    double m_alpha_1;
    double m_alpha_2;
    double m_beta_1;
    double m_beta_2;
    double m_X_0;
    double m_X_1;
    double m_Y_0;
    double m_Y_1;
    double m_Z_0;
    double m_Z_1;

    arma::mat m_x;
    arma::vec m_js;
    std::vector<vector<int>> m_gridNeighbours;
    std::vector<vector<int>> m_particlesInGridPoint;
    std::vector<double> m_gridPointPorosity;

    unsigned seed;
    std::default_random_engine m_generator;
    std::uniform_real_distribution<double> m_distribution_x;
    std::uniform_real_distribution<double> m_distribution_y;
    std::uniform_real_distribution<double> m_distribution_z;


    double m_dx;
    double m_dy;
    double m_dz;
    double m_DX;
    double m_DY;
    double m_DZ;

    bool m_periodic_x;
    bool m_periodic_y;
    bool m_periodic_z;

    double m_optimalGridSpacing;

    // Domain variables
    int m_dim = 3;
    double m_nx, m_ny, m_nz;
    double m_gridSpacing_x, m_gridSpacing_y, m_gridSpacing_z;

    int m_imageResolution;
    int m_testSaveFreq;

    int findGridId(const arma::vec3 & r_i);
    void checkBoundaries();

    int m_openmp_threads;
//    int m_globalGridSize = 10.01;
//    int m_globalGridSize = 6.01;
//    int m_globalGridSize = 2.01;
    int m_globalGridSize = 1.51;
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
