#include "meshgenerator.h"

//------------------------------------------------------------------------------
mg::MeshGenerator::MeshGenerator()
{

}
//------------------------------------------------------------------------------
mg::MeshGenerator::MeshGenerator(mg::Parameters parameters):
    param(parameters)
{

    if(parameters.dataType == "binvox") {
        readBinvox();
    } else if(parameters.dataType == "images") {
            readImages();
    } else if(parameters.dataType == "image") {
        readImage();
    } else {
        std::cerr << parameters.dataType << " is not a valid datatype" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    m_nParticles = parameters.nParticles;
    m_samples = parameters.q;
    m_steps = parameters.threshold;

    m_alpha_1 = parameters.alpha_1;
    m_alpha_2 = parameters.alpha_2;
    m_beta_1 = parameters.beta_1;
    m_beta_2 = parameters.beta_2;
    m_X_0 = 0;
    m_X_1 = (double)m_w/m_w;
    m_Y_0 = 0;
    m_Y_1 = (double)m_h/m_w;
    m_Z_0 = 0;
    m_Z_1 = parameters.z_scale*(double)(m_d)/m_w;

    m_x = arma::randu(3,m_nParticles);
    m_js = arma::ones(m_nParticles);

    seed = std::chrono::system_clock::now().time_since_epoch().count();
    cout << "Seed:" << seed << endl;
    m_generator = std::default_random_engine();
    m_generator.seed(seed);

    m_distribution_x = std::uniform_real_distribution<double> (m_X_0, m_X_1);
    m_distribution_y = std::uniform_real_distribution<double> (m_Y_0, m_Y_1);
    m_distribution_z = std::uniform_real_distribution<double> (m_Z_0, m_Z_1);
//    cout << "Test generation " <<  m_distribution_x(m_generator) << " " <<  m_distribution_y(m_generator) << " " <<  m_distribution_z(m_generator) << endl;
    m_dx  = (m_X_1 - m_X_0)/m_w;
    m_dy  = (m_Y_1 - m_Y_0)/m_h;
    m_dz  = (m_Z_1 - m_Z_0)/m_d;
    m_DX = (m_X_1 - m_X_0);
    m_DY = (m_Y_1 - m_Y_0);
    m_DZ = (m_Z_1 - m_Z_0);

    m_periodic_x = parameters.periodic_x;
    m_periodic_y = parameters.periodic_y;
    m_periodic_z = parameters.periodic_z;

    m_imageResolution = parameters.imageResolution;
    m_testSaveFreq = parameters.testSaveFreq;

    setDomainSize(m_globalGridSize, m_solidFraction);

    m_openmp_threads = parameters.openmp_threads;
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::initializeFromImage()
{
    const int max_number_of_tries = 60000;
    const double estimatedParticleLength = 1.3*pow(m_DX*m_DY*m_DZ*m_solidFraction/m_nParticles, 1./3.);
    Grid grid(m_X_0, m_X_1, m_Y_0, m_Y_1, m_Z_0, m_Z_1, estimatedParticleLength);
    grid.createGrid();

    computePorosity(grid, m_img_data);
    const vector<double> & porosities = grid.getGridProperties("porosity");

    const vector<arma::vec3> & gridCenters= grid.getCenters();

    //--------------------------------------
    const vector<double> & gridSpacing = grid.getGridspacing();
    const int nElements = grid.nElements();
    double sumParticles = 0;
    int nPlaced = 0;

    std::uniform_real_distribution<double> distribution_01 =  std::uniform_real_distribution<double> (0., 1.);
    double r_i[3] = {0,0,0};

//    cout << "Test generation2 " <<  m_distribution_x(m_generator) << " " <<  m_distribution_y(m_generator) << " " <<  m_distribution_z(m_generator) << endl;
    cout << "placing particles:" << endl;
    for(int gId=0; gId<nElements; gId++) {
        const double estimated_nParticles = porosities[gId]*m_nParticles;
        sumParticles += estimated_nParticles;

        int nParticles = floor(estimated_nParticles);
        double prop = 0.99*(estimated_nParticles-nParticles);
        if(prop < 0.05)
            prop = 0;

        if(distribution_01(m_generator) < prop)
            nParticles++;

        bool notFound = false;
        if(nParticles > 0) {
            const arma::vec3 center = gridCenters[gId];
            double x0 = center[0] - 0.5*gridSpacing[0];
            double x1 = center[0] + 0.5*gridSpacing[0];
            double y0 = center[1] - 0.5*gridSpacing[1];
            double y1 = center[1] + 0.5*gridSpacing[1];
            double z0 = center[2] - 0.5*gridSpacing[2];
            double z1 = center[2] + 0.5*gridSpacing[2];

            m_distribution_x = std::uniform_real_distribution<double> (x0, x1);
            m_distribution_y = std::uniform_real_distribution<double> (y0, y1);
            m_distribution_z = std::uniform_real_distribution<double> (z0, z1);

            for(int i=nPlaced; i<nPlaced+nParticles; i++) {
                int counter = 0;
                do {
                    r_i[0] = m_distribution_x(m_generator);
                    r_i[1] = m_distribution_y(m_generator);
                    r_i[2] = m_distribution_z(m_generator);
                    counter++;
                    if(counter>max_number_of_tries) {
                        notFound = true;
                        break;
                    }
                } while (m_img_data(r_i[0]/m_dx, r_i[1]/m_dy, r_i[2]/m_dz) > 0);
                m_x(0,i) = r_i[0];
                m_x(1,i) = r_i[1];
                m_x(2,i) = r_i[2];
            }
            if(notFound) {
                cerr << "No particles places in section " << gId << endl;
            } else {
                nPlaced += nParticles;
            }
        }
        if(nPlaced%200 == 0)
            printProgress(double(nPlaced)/m_nParticles);
    }

    cout << "sumParticles:" << sumParticles << " nPlaced:" << nPlaced  << endl;

    //--------------------------------------

    m_distribution_x = std::uniform_real_distribution<double> (m_X_0, m_X_1);
    m_distribution_y = std::uniform_real_distribution<double> (m_Y_0, m_Y_1);
    m_distribution_z = std::uniform_real_distribution<double> (m_Z_0, m_Z_1);

    // Randomly trying points within the image.
    for(int i=nPlaced; i<m_nParticles; i++) {
        do {
            r_i[0] = m_distribution_x(m_generator);
            r_i[1] = m_distribution_y(m_generator);
            r_i[2] = m_distribution_z(m_generator);
        } while (m_img_data(r_i[0]/m_dx, r_i[1]/m_dy, r_i[2]/m_dz) > 0);
        m_x(0, i) = r_i[0];
        m_x(1, i) = r_i[1];
        m_x(2, i) = r_i[2];
    }

    std::cout << "Initialization from image complete" << std::endl;
}
//------------------------------------------------------------------------------
arma::mat mg::MeshGenerator::createMesh()
{
    cout <<  "createDomainGrid()" << endl;
    createDomainGrid();
    cout <<  "initializeFromImage()" << endl;
    initializeFromImage();
    cout <<  "mapParticlesToGrid()" << endl;
    mapParticlesToGrid();
    cout <<  "starting" << endl;

//    computeGridPorosity();

    std::uniform_real_distribution<double> distribution_rand_particle(m_X_0, m_nParticles);
    // Sampling the image Monte Carlo style and adjusting the point centers
    // untill convergence.


    vector<array<double, 4>> neighbours(m_nParticles);
#ifdef FORCE_OMP_CPU
        omp_set_num_threads(m_openmp_threads);
#endif
#pragma omp parallel for
    for(int i=0; i<m_nParticles; i++) {
        array<double, 4> &dun = neighbours[i];
        dun[0] = 0;
        dun[1] = 0;
        dun[2] = 0;
        dun[3] = 0;
    }

    for (int k=0; k<m_steps;k++) {
        printProgress(double(k)/m_steps);

        checkBoundaries();
        mapParticlesToGrid();

        if(param.testingSave  && k % m_testSaveFreq == 0) {
            save_image_and_xyz(param.savePath + "/alg1", k);
            calculateRadialDistribution(k);
        }

        if(k % param.redistributionFrequency == 0) {
            // Picking nRandom points for redistribution
            for(int iterations=0; iterations < param.nRedistributedPoints; iterations++) {
                int random_particle = distribution_rand_particle(m_generator);
                do {
                    m_x(0, random_particle) = m_distribution_x(m_generator);
                    m_x(1, random_particle) = m_distribution_y(m_generator);
                    m_x(2, random_particle) = m_distribution_z(m_generator);
            }while(m_img_data(m_x(0, random_particle)/m_dx, m_x(1, random_particle)/m_dy, m_x(2, random_particle)/m_dz) > 0);
            }
        }

//        unordered_map <int, vector <arma::vec3>> neighbours;

#ifdef FORCE_OMP_CPU
        omp_set_num_threads(m_openmp_threads);
#endif
#pragma omp parallel for
        for(int r=0; r<m_samples; r++) {
            arma::vec3 y_r;
            double maxLen = numeric_limits<double>::max();
            int indexMax = -1;

            do {
                y_r[0] = m_distribution_x(m_generator);
                y_r[1] = m_distribution_y(m_generator);
                y_r[2] = m_distribution_z(m_generator);
            } while(m_img_data(y_r(0)/m_dx, y_r(1)/m_dy, y_r(2)/m_dz) > 0);

            arma::vec3 y_tmp = y_r;
            const int gId = findGridId(y_tmp);
            // Finding the closest voronoi center
            //------------------------------------------------------------------
            // Checking this gridpoint
            //------------------------------------------------------------------

            for(int k:m_particlesInGridPoint[gId]) {
                arma::vec3 y_r_copy = y_r;
                arma::vec3 x_k = y_r - m_x.col(k);

                if(m_periodic_x) {
                    if(x_k(0) > 0.5*m_DX) {
                        x_k(0) -= m_DX;
                        y_r_copy(0) -= m_DX;
                    }else if(x_k(0) < -0.5*m_DX) {
                        x_k(0) += m_DX;
                        y_r_copy(0) += m_DX;
                    }
                }

                if(m_periodic_y) {
                    if(x_k(1) > 0.5*m_DY) {
                        x_k(1) -= m_DY;
                        y_r_copy(1) -= m_DY;
                    }else if(x_k(1) < -0.5*m_DY) {
                        x_k(1) += m_DY;
                        y_r_copy(1) += m_DY;
                    }
                }

                if(m_periodic_z) {
                    if(x_k(2) > 0.5*m_DZ) {
                        x_k(2) -= m_DZ;
                        y_r_copy(2) -= m_DZ;
                    }else if(x_k(2) < -0.5*m_DZ) {
                        x_k(2) += m_DZ;
                        y_r_copy(2) += m_DZ;
                    }
                }
#if NORM_N==2
                 double dr_rk = x_k(0)*x_k(0) + x_k(1)*x_k(1) + x_k(2)*x_k(2);
#else
                const double dr_rk = pow(x_k(0), NORM_N) + pow(x_k(1), NORM_N) + pow(x_k(2), NORM_N);
#endif
                if(dr_rk < maxLen) {
                    maxLen = dr_rk;
                    indexMax = k;
                    y_tmp = y_r_copy;
                }
            }

            //------------------------------------------------------------------
            // Checking neighbouring gridpoint
            //------------------------------------------------------------------
            for(int gridNeighbour:m_gridNeighbours[gId]) {
                for(int k:m_particlesInGridPoint[gridNeighbour]) {
                    arma::vec3 y_r_copy = y_r;
                    arma::vec3 x_k = y_r - m_x.col(k);

                    if(m_periodic_x) {
                        if(x_k(0) > 0.5*m_DX) {
                            x_k(0) -= m_DX;
                            y_r_copy(0) -= m_DX;
                        }else if(x_k(0) < -0.5*m_DX) {
                            x_k(0) += m_DX;
                            y_r_copy(0) += m_DX;
                        }
                    }

                    if(m_periodic_y) {
                        if(x_k(1) > 0.5*m_DY) {
                            x_k(1) -= m_DY;
                            y_r_copy(1) -= m_DY;
                        }else if(x_k(1) < -0.5*m_DY) {
                            x_k(1) += m_DY;
                            y_r_copy(1) += m_DY;
                        }
                    }

                    if(m_periodic_z) {
                        if(x_k(2) > 0.5*m_DZ) {
                            x_k(2) -= m_DZ;
                            y_r_copy(2) -= m_DZ;
                        }else if(x_k(2) < -0.5*m_DZ) {
                            x_k(2) += m_DZ;
                            y_r_copy(2) += m_DZ;
                        }
                    }

#if NORM_N==2
                    double dr_rk = x_k(0)*x_k(0) + x_k(1)*x_k(1) + x_k(2)*x_k(2);
#else
                    const double dr_rk = pow(x_k(0), NORM_N) + pow(x_k(1), NORM_N) + pow(x_k(2), NORM_N);
#endif
                    if(dr_rk < maxLen) {
                        maxLen = dr_rk;
                        indexMax = k;
                        y_tmp = y_r_copy;
                    }
                }
            }
            //------------------------------------------------------------------
            // Storing the result
#pragma omp critical
            {
                if(indexMax >= 0) {
                    array<double, 4> &du = neighbours[indexMax];
                    du[0] += y_tmp[0];
                    du[1] += y_tmp[1];
                    du[2] += y_tmp[2];
                    du[3] += 1;
                }
            }
        }


        arma::vec3 u_r;
//#ifdef FORCE_OMP_CPU
//        omp_set_num_threads(openmp_threads);
//#pragma omp parallel for private(u_r)
        for(int i=0; i<m_nParticles; i++) {
            array<double, 4> &dun = neighbours[i];

            if(dun[3] <= 0)
                continue;
            double j = m_js(i);
            const arma::vec3 &x_i =  m_x.col(i);
            u_r[0] = dun[0];
            u_r[1] = dun[1];
            u_r[2] = dun[2];

            u_r /= dun[3];

            m_x.col(i) = ((m_alpha_1*j + m_beta_1)*x_i + (m_alpha_2*j + m_beta_2)*u_r)/(j+1);
            m_js(i) += 1;

            dun[0] = 0;
            dun[1] = 0;
            dun[2] = 0;
            dun[3] = 0;

            /*
        const vector <arma::vec3> neigh = neighbours[i];

        if(neigh.size() > 0) {
            const arma::vec3 x_i =  x.col(i);
            double j = js(i);
            // Finding the average of the random points
            arma::vec3 u_r = arma::zeros(3);
            for (uint l=0; l<neigh.size(); l++)
            {
                const arma::vec3& yr = neigh[l];
                u_r += yr;
            }
            u_r /= neigh.size();

            x.col(i) = ((alpha_1*j + beta_1)*x_i + (alpha_2*j + beta_2)*u_r)/(j+1);
            js(i) += 1;
        }*/
        }
    }

//    computeGridPorosity();
    return m_x;
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::createDomainGrid()
{
    std::vector<int> pluss_minus = {-1, 0, 1};
    m_gridNeighbours = std::vector<vector<int>> (m_nx*m_ny*m_nz, std::vector<int>(0));
    m_particlesInGridPoint.clear();
    m_particlesInGridPoint = std::vector<vector<int>>(m_nx*m_ny*m_nz, std::vector<int>(0));

    for(int x_=0;x_<m_nx;x_++) {
        for(int y_=0;y_<m_ny;y_++) {
            for(int z_=0;z_<m_nz;z_++) {
                int id_ = z_ + m_nz*y_ + m_nz*m_ny*x_;

                for(int id_x:pluss_minus) {
                    id_x += x_;
                    if(id_x < 0) {
                        if(m_periodic_x)
                            id_x = m_nx-1;
                        else
                            continue;
                    }
                    if(id_x >= m_nx) {
                        if(m_periodic_x)
                            id_x = 0;
                        else
                            continue;
                    }

                    for(int id_y:pluss_minus) {
                        id_y += y_;
                        if(id_y < 0) {
                            if(m_periodic_y)
                                id_y = m_ny-1;
                            else
                                continue;
                        }
                        if(id_y >= m_ny) {
                            if(m_periodic_y)
                                id_y = 0;
                            else
                                continue;
                        }

                        for(int id_z:pluss_minus) {
                            id_z += z_;
                            if(id_z < 0) {
                                if(m_periodic_z)
                                    id_z = m_nz-1;
                                else
                                    continue;
                            }
                            if(id_z >= m_nz) {
                                if(m_periodic_z)
                                    id_z = 0;
                                else
                                    continue;
                            }

                            int id_neighbour = id_z + m_nz*id_y + m_nz*m_ny*id_x;
                            if(id_neighbour == id_)
                                continue;
                            m_gridNeighbours[id_].push_back(id_neighbour);
                        }
                    }
                }
            }
        }
    }
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::mapParticlesToGrid()
{
    // Emptying the grid
    for(vector<int> &t:m_particlesInGridPoint) {
        t.clear();
    }
    // Placing all particles in the grid

#ifdef FORCE_OMP_CPU
        omp_set_num_threads(m_openmp_threads);
#endif
#pragma omp parallel for
    for(int i=0; i<m_nParticles; i++) {
        const arma::vec3 & r_i = m_x.col(i);
#pragma omp critical
        m_particlesInGridPoint[findGridId(r_i)].push_back(i);
    }
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::save_image_and_xyz(string base, int nr)
{
    // Bounds check
    checkBoundaries();
    mapParticlesToGrid();

    arma::vec areas = arma::zeros(m_nParticles);

    int resolution_x = m_X_1*m_imageResolution;
    int resolution_y = m_Y_1*m_imageResolution;
    int resolution_z = m_Z_1*m_imageResolution;
    int pix_hole = 0;
    //--------------------------------------------------------------------------
    // Creating a Voronoi image and computing the areas
    //--------------------------------------------------------------------------


#ifdef FORCE_OMP_CPU
        omp_set_num_threads(m_openmp_threads);
#endif
#pragma omp parallel for
    for (int i=0; i<resolution_x;i++) {
        arma::vec3 r_img;
        arma::vec3 x_k;

        for (int j=0; j<resolution_y;j++) {
            r_img[0] = m_X_1*i/(resolution_x);
            r_img[1] = m_Y_1*j/(resolution_y);

            for (int k=0; k<resolution_z;k++) {
                r_img[2] = m_Z_1*k/(resolution_z);

                if(m_img_data(r_img(0)/m_dx, r_img(1)/m_dy, r_img(2)/m_dz) > 0) {
//                    if(m_img_data(r_img(1)/m_dy, r_img(0)/m_dx, r_img(2)/m_dz) > 0) {
                    pix_hole++;
                    continue;
                }
                double maxLen = numeric_limits<double>::max();
                int indexMax = -1;
                int gId = findGridId(r_img);

                // Finding the closest voronoi center
                //------------------------------------------------------------------
                // Checking this gridpoint
                //------------------------------------------------------------------
                for(int k:m_particlesInGridPoint[gId]) {
                    x_k = r_img - m_x.col(k);

                    if(m_periodic_x) {
                        if(x_k(0) > 0.5*m_DX) {
                            x_k(0) -= m_DX;
                        }else if(x_k(0) < -0.5*m_DX) {
                            x_k(0) += m_DX;
                        }
                    }

                    if(m_periodic_y) {
                        if(x_k(1) > 0.5*m_DY) {
                            x_k(1) -= m_DY;
                        }else if(x_k(1) < -0.5*m_DY) {
                            x_k(1) += m_DY;
                        }
                    }

                    if(m_periodic_z) {
                        if(x_k(2) > 0.5*m_DZ) {
                            x_k(2) -= m_DZ;
                        }else if(x_k(2) < -0.5*m_DZ) {
                            x_k(2) += m_DZ;
                        }
                    }
#if NORM_N==2
                    double dr_rk = (x_k(0)*x_k(0) + x_k(1)*x_k(1)) + x_k(2)*x_k(2);
#else
                    const double dr_rk = pow(x_k(0), NORM_N) + pow(x_k(1), NORM_N) + pow(x_k(2), NORM_N);
#endif
                    if(dr_rk < maxLen) {
                        maxLen = dr_rk;
                        indexMax = k;
                    }
                }

                //------------------------------------------------------------------
                // Checking neighbouring gridpoint
                //------------------------------------------------------------------
                for(int gridNeighbour:m_gridNeighbours[gId]) {
                    for(int k:m_particlesInGridPoint[gridNeighbour]) {
                        x_k = r_img - m_x.col(k);

                        if(m_periodic_x) {
                            if(x_k(0) > 0.5*m_DX) {
                                x_k(0) -= m_DX;
                            }else if(x_k(0) < -0.5*m_DX) {
                                x_k(0) += m_DX;
                            }
                        }

                        if(m_periodic_y) {
                            if(x_k(1) > 0.5*m_DY) {
                                x_k(1) -= m_DY;
                            }else if(x_k(1) < -0.5*m_DY){
                                x_k(1) += m_DY;
                            }
                        }

                        if(m_periodic_z) {
                            if(x_k(2) > 0.5*m_DZ) {
                                x_k(2) -= m_DZ;
                            }else if(x_k(2) < -0.5*m_DZ){
                                x_k(2) += m_DZ;
                            }
                        }

#if NORM_N==2
                    double dr_rk = (x_k(0)*x_k(0) + x_k(1)*x_k(1)) + x_k(2)*x_k(2);
#else
                    const double dr_rk = pow(x_k(0), NORM_N) + pow(x_k(1), NORM_N) + pow(x_k(2), NORM_N);
#endif
                        if(dr_rk < maxLen) {
                            maxLen = dr_rk;
                            indexMax = k;
                        }
                    }
                }

                if(indexMax != -1) {
                    areas[indexMax] += 1.0;
                }
            }
        }
    }


    //--------------------------------------------------------------------------
    // Saving xyz-file with volume
    //--------------------------------------------------------------------------
    string fileName;
    if(nr == -1)
        fileName =  base + ".xyz";
    else
        fileName =  base + "_" + to_string(nr) + ".xyz";
    ofstream outStream(fileName.c_str());

    outStream << m_x.n_cols << endl;
    outStream << "# id x y z volume" << endl;
    double dxdydz = (m_X_1 - m_X_0)*(m_Y_1 - m_Y_0)*(m_Z_1 - m_Z_0);
    double s_volume = 0;
//    double total_pix = resolution_x*resolution_y*resolution_z - pix_hole;
    double total_pix = resolution_x*resolution_y*resolution_z;
//    double optimalPackingOfSpheres = 0.68;
    double optimalPackingOfSpheres = 1.0;

    for (int i=0; i<m_nParticles;i++) {
        double volume = optimalPackingOfSpheres*dxdydz*areas[i]/(total_pix);
        const arma::vec3& r = m_x.col(i);
        outStream << i << "\t" << r(0) << "\t" << r(1) << "\t" << r(2) << "\t" << volume << std::endl;
        s_volume += volume;
    }
    outStream.close();
    cout << fileName << endl;

    cout << "Volumefraction:" << s_volume/dxdydz << endl;
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::setDomainSize(double spacing, double solidFraction)
{
    // Setting the grid size

    double estimatedParticleLength = pow(m_DX*m_DY*m_DZ*solidFraction/m_nParticles, 1./3.);
    double gridSpacing = spacing*estimatedParticleLength;

    m_nx = floor((m_DX)/gridSpacing);
    m_ny = floor((m_DY)/gridSpacing);
    m_nz = floor((m_DZ)/gridSpacing);

    if(m_nx == 0)
        m_nx = 1;
    if(m_ny == 0)
        m_ny = 1;
    if(m_nz == 0)
        m_nz = 1;

    m_gridSpacing_x = m_DX/m_nx;
    m_gridSpacing_y = m_DY/m_ny;
    m_gridSpacing_z = m_DZ/m_nz;
}
//------------------------------------------------------------------------------
double mg::MeshGenerator::calculateRadialDistribution(int nr)
{
    std::cout << "Calculating histogram" << std::endl;
//    setDomainSize(7.5*m_globalGridSize, m_solidFraction);
//    createDomainGrid();
    checkBoundaries();
    std::cout << "pre mapParticlesToGrid" << std::endl;
    mapParticlesToGrid();
    std::cout << "Post mapParticlesToGrid" << std::endl;

    int nBins = 400;
    double maxLength = 2.0*m_gridSpacing_x;
    double histSpacing = maxLength/nBins;
    vector<int> histogram(nBins, 0);


#ifdef FORCE_OMP_CPU
        omp_set_num_threads(m_openmp_threads);
#endif
#pragma omp parallel for
    for(int i=0; i<m_nParticles; i++)
    {
        const arma::vec3 & r_i = m_x.col(i);
        int gId = findGridId(r_i);

        //------------------------------------------------------------------
        // Checking this gridpoint
        //------------------------------------------------------------------

        for(int j:m_particlesInGridPoint[gId]) {
            if(j == i)
                continue;
            arma::vec3 r_ij = r_i - m_x.col(j);

            if(m_periodic_x) {
                if(r_ij(0) > 0.5*m_DX) {
                    r_ij(0) -= m_DX;
                }else if(r_ij(0) < -0.5*m_DX) {
                    r_ij(0) += m_DX;
                }
            }

            if(m_periodic_y) {
                if(r_ij(1) > 0.5*m_DY){
                    r_ij(1) -= m_DY;
                }else if(r_ij(1) < -0.5*m_DY){
                    r_ij(1) += m_DY;
                }
            }

            if(m_periodic_z) {
                if(r_ij(2) > 0.5*m_DZ) {
                    r_ij(2) -= m_DZ;
                }else if(r_ij(2) < -0.5*m_DZ) {
                    r_ij(2) += m_DZ;
                }
            }

#if NORM_N==2
            double dr_rk = sqrt(r_ij(0)*r_ij(0) + r_ij(1)*r_ij(1) + r_ij(2)*r_ij(2));
#else
            const double dr_rk = pow(pow(r_ij(0), NORM_N) + pow(r_ij(1), NORM_N) + pow(r_ij(2), NORM_N), 1./NORM_N);
#endif
            if(dr_rk > maxLength)
                continue;
            int id = dr_rk/histSpacing;
            histogram[id]++;
        }

        //------------------------------------------------------------------
        // Checking neighbouring gridpoint
        //------------------------------------------------------------------
        for(int gridNeighbour:m_gridNeighbours[gId]) {
            for(int k:m_particlesInGridPoint[gridNeighbour]) {
                arma::vec3 r_ij = r_i - m_x.col(k);

                if(m_periodic_x) {
                    if(r_ij(0) > 0.5*m_DX) {
                        r_ij(0) -= m_DX;
                    }else if(r_ij(0) < -0.5*m_DX) {
                        r_ij(0) += m_DX;
                    }
                }

                if(m_periodic_y) {
                    if(r_ij(1) > 0.5*m_DY) {
                        r_ij(1) -= m_DY;
                    }else if(r_ij(1) < -0.5*m_DY) {
                        r_ij(1) += m_DY;
                    }
                }

                if(m_periodic_z) {
                    if(r_ij(2) > 0.5*m_DZ) {
                        r_ij(2) -= m_DZ;
                    }else if(r_ij(2) < -0.5*m_DZ) {
                        r_ij(2) += m_DZ;
                    }
                }
#if NORM_N==2
            double dr_rk = sqrt(r_ij(0)*r_ij(0) + r_ij(1)*r_ij(1) + r_ij(2)*r_ij(2));
#else
            const double dr_rk = pow(pow(r_ij(0), NORM_N) + pow(r_ij(1), NORM_N) + pow(r_ij(2), NORM_N), 1./NORM_N);
#endif
                if(dr_rk > maxLength)
                    continue;
                int id = dr_rk/histSpacing;
                histogram[id]++;
            }
        }
        //------------------------------------------------------------------
    }
    string fileName;
    if(nr == -1)
        fileName = param.savePath + "/histogram.hist";
    else
        fileName = param.savePath + "/histogram_" + to_string(nr) + ".hist";
    ofstream outStream(fileName.c_str());

    // Finding the optimal spacing between the particles
    int maxIndex = -1;
    double maxValue = -999;
    m_optimalGridSpacing = 0;

    for(int i=0;i<nBins; i++)
    {
        double r = i*histSpacing;
        double r2 = (i + 1.0)*histSpacing;
        double val = histogram[i]/((pow(r2, 3) - pow(r,3)));
        if(val > maxValue)
        {
            maxIndex = i;
            maxValue = val;
        }
        outStream << r + 0.5*histSpacing << "\t" << val << std::endl;
    }
    outStream.close();

    m_optimalGridSpacing = (maxIndex + 0.5)*histSpacing;

    // Resetting the grid
    setDomainSize(m_globalGridSize, m_solidFraction);
    createDomainGrid();
    mapParticlesToGrid();

    return m_optimalGridSpacing;
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::writeConfiguration()
{
    std::cout << "Writing configuration" << std::endl;
    string fileName = param.savePath + "/configuration.cfg";
    ofstream outStream(fileName.c_str());
    outStream.setf(ios::scientific);
    outStream.precision(5);

    outStream << "nParticles = " << m_nParticles << std::endl;
    outStream << "spacing = " << m_optimalGridSpacing << std::endl;
    int n_x = floor((m_X_1 - m_X_0)/m_optimalGridSpacing);
    int n_y = floor((m_Y_1 - m_Y_0)/m_optimalGridSpacing);
    int n_z = floor((m_Z_1 - m_Z_0)/m_optimalGridSpacing);
    outStream << "latticePoints = [" << n_x<< ", " << n_y << ", " << n_z << "]"
              << std::endl;
    outStream << "boundaries = [" << m_X_0 << ", " << m_X_1 <<  ", "
              << m_Y_0 << ", " << m_Y_1 << ", "
              << m_Z_0 << ", " << m_Z_1 << "]"
              << std::endl;
    outStream << "dxdydx = [" << m_dx << ", "
              << m_dy << ", "
              << m_dz << "]"
              << std::endl;

    outStream << "periodic = [";
    if(m_periodic_x)
        outStream << "1, ";
    else
        outStream << "0, ";

    if(m_periodic_y)
        outStream << "1, ";
    else
        outStream << "0, ";

    if(m_periodic_z)
        outStream << "1";
    else
        outStream << "0";

    outStream << "]" << std::endl;

    outStream.close();
}
//------------------------------------------------------------------------------
int mg::MeshGenerator::findGridId(const arma::vec3 &r)
{
    int id_x = (r(0) - m_X_0)/m_gridSpacing_x;
    int id_y = (r(1) - m_Y_0)/m_gridSpacing_y;
    int id_z = (r(2) - m_Z_0)/m_gridSpacing_z;

    // Boundary checks
    if(id_x > m_nx)
        id_x = m_nx - 1;
    else if(id_x < 0)
        id_x = 0;

    if(id_y > m_ny)
        id_y = m_ny - 1;
    else if(id_y < 0)
        id_y = 0;

    if(id_z > m_nz)
        id_z = m_nz - 1;
    else if(id_z < 0)
        id_z = 0;

    return id_z + m_nz*id_y + m_ny*m_nz*id_x;
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::checkBoundaries()
{
    if(m_periodic_x) {
#pragma omp parallel for
        for(int k=0;k<m_nParticles; k++) {
            if(m_x(0,k) < m_X_0)
                m_x(0,k) += m_DX;
            if(m_x(0,k) >= m_X_1)
                m_x(0,k) -= m_DX;
        }
    }

    if(m_periodic_y) {
#pragma omp parallel for
        for(int k=0;k<m_nParticles; k++)
        {
            if(m_x(1,k) < m_Y_0)
                m_x(1,k) += m_DY;
            if(m_x(1,k) >= m_Y_1)
                m_x(1,k) -= m_DY;
        }
    }

    if(m_periodic_z) {
#pragma omp parallel for
        for(int k=0;k<m_nParticles; k++)
        {
            if(m_x(2,k) < m_Z_0)
                m_x(2,k) += m_DZ;
            if(m_x(2,k) >= m_Z_1)
                m_x(2,k) -= m_DZ;
        }
    }
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::readImages()
{
    int zeros = param.leading_zeros;
    m_d = param.imgRange[1] - param.imgRange[0] + 1;

    // First case
    string fileName = param.dataPath + "/" + param._baseName;
    string n_string = to_string(param.imgRange[0]);
    fileName += string(zeros - n_string.length(), '0') + n_string;
    fileName += "." + param.imageFormat;
    CImg<double> f_image = CImg<double> (fileName.c_str());

    m_h = f_image.height();
    m_w = f_image.width();
    m_img_data = cube_type(m_h, m_w, m_d);
    const double max_value = 255.;
    double solidFraction = 0;

    // The rest of the images
    int k = 0;
    for(int z_range = param.imgRange[1]; z_range >= param.imgRange[0]; z_range-- ) {
        fileName = param.dataPath + "/" + param._baseName;
        n_string = to_string(z_range);
        fileName += string(zeros - n_string.length(), '0') + n_string;
        fileName += "." + param.imageFormat;
        CImg<double> image = CImg<double> (fileName.c_str());

        for(int j=0;j<m_w;j++) {
            for(int i=0;i<m_h;i++) {
                const double weight = 1. - image(j, i, 0, 0)/(max_value);
                m_img_data(i, j, k) = weight;
                if(weight == 0) {
                    solidFraction += 1;
                }
            }
        }
        k++;
    }
    m_solidFraction = solidFraction/(m_w*m_h*m_d);

    std::cout << "Images series read from file" << std::endl;
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::readImage()
{
    CImg<double> image = CImg<double> (param.dataPath.c_str());

    m_h = image.height();
    m_w = image.width();
    m_d = param.z_depth;

    m_img_data = cube_type(m_h, m_w, m_d);
    for(int j=0;j<m_w;j++) {
        for(int i=0;i<m_h;i++) {
            for(int k=0;k<m_d;k++) {
                m_img_data(i,j, k) = 1. - image(j,i,0,0)/(255.0);
            }
        }
    }
    std::cout << "Image read from file" << std::endl;
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::readBinvox()
{
    // From http://www.cs.princeton.edu/~min/binvox/read_binvox.html

    typedef unsigned char byte;
    std::ifstream *input = new std::ifstream(param.dataPath.c_str(), std::ios::in | std::ios::binary);

    //-------------------
    // binvox header
    //-------------------
    string line;
    *input >> line;  // #binvox
    if (line.compare("#binvox") != 0) {
      cout << "Error: first line reads [" << line << "] instead of [#binvox]" << endl;
      delete input;
      std::exit(EXIT_FAILURE);
    }

    int version;
    *input >> version;

    int depth, height, width;
    depth = -1;
    int done = 0;

    while(input->good() && !done) {
        *input >> line;
        if (line.compare("data") == 0) done = 1;
        else if (line.compare("dim") == 0) {
            *input >> depth >> height >> width;
        } else {
            cout << "  unrecognized keyword [" << line << "], skipping" << endl;
            char c;
            do
            {  // skip until end of line
                c = input->get();
            } while(input->good() && (c != '\n'));

        }
    }

    if (!done) {
      cout << "  error reading header" << endl;
      std::exit(EXIT_FAILURE);
    }

    if (depth == -1) {
      cout << "  missing dimensions in header" << endl;
      std::exit(EXIT_FAILURE);
    }

    int size = width * height * depth;
    byte *voxels = new byte[size];
    if (!voxels) {
      cout << "  error allocating memory" << endl;
      std::exit(EXIT_FAILURE);
    }

    //-------------------
    // read voxel data
    //-------------------
    byte value;
    byte count;
    int index = 0;
    int end_index = 0;
    int nr_voxels = 0;

    input->unsetf(ios::skipws);  // need to read every byte now (!)
    *input >> value;  // read the linefeed char

    while((end_index < size) && input->good()) {
      *input >> value >> count;

      if (input->good()) {
        end_index = index + count;

        if (end_index > size)
            std::exit(EXIT_FAILURE);

        for(int i=index; i < end_index; i++) {
            voxels[i] = value;
        }
        if (value)
            nr_voxels += count;
        index = end_index;
      }  // if file still ok

    }  // while

    input->close();

    double wxh = width * height;

    int x_min = 2*height;
    int x_max = 0;
    int y_min = 2*height;
    int y_max = 0;
    int z_min = 2*height;
    int z_max = 0;
    for(int j=0;j<width;j++)
    {
        for(int i=0;i<height;i++)
        {
            for(int k=0;k<depth;k++)
            {
                int index = j * wxh + i * width + k;  // wxh = width * height = d * d
                if(voxels[index])
                {
                    if(i<x_min)
                        x_min = i;
                    if(i>x_max)
                        x_max = i;

                    if(j<y_min)
                        y_min = j;
                    if(j>y_max)
                        y_max = j;

                    if(k<z_min)
                        z_min = k;
                    if(k>z_max)
                        z_max = k;
                }
            }
        }
    }

    m_h = x_max - x_min;
    m_w = y_max - y_min;
    m_d = z_max - z_min;
    m_img_data = cube_type(m_h, m_w, m_d);

    for(int j=y_min; j<y_max; j++) {
        for(int i=x_min;i<x_max;i++) {
            for(int k=z_min;k<z_max;k++) {
                int index = j * wxh + i * width + k;
                if(voxels[index])
                    m_img_data(i, j, k) = 1;
//                if(!voxels[index])
//                img_data(i, j, k) = 1;
            }
        }
    }

    cout << "binvox file read" << endl;
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::printProgress(double progress)
{
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int j = 0; j < barWidth; ++j) {
        if (j < pos)
            std::cout << "=";
        else if (j == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}
//------------------------------------------------------------------------------
void mg::MeshGenerator::computeGridPorosity()
{
    const int nElements = m_nx*m_ny*m_nz;
    m_gridPointPorosity = std::vector<double>(nElements, 0);
    std::vector<int> nPointsInGridPoint = std::vector<int>(nElements, 0);
    arma::vec3 r = {0,0,0};
    double solidfraction = 0;

    for(int i=0;i<m_w;i++) {
        for(int j=0;j<m_h;j++) {
            for(int k=0;k<m_d;k++) {
                r = {j*m_dx, i*m_dy, k*m_dz};
                const int gId = findGridId(r);
                nPointsInGridPoint[gId] += 1;

                if (m_img_data(i, j, k) == 0) {
                    m_gridPointPorosity[gId] += 1;
                    solidfraction  += 1;
                }
            }
        }
    }
    solidfraction  /= (m_w*m_h*m_d);
    const double particleVolume = m_DX*m_DZ*m_DZ/(double)m_nParticles;
    const double volumegridElement = m_gridSpacing_x*m_gridSpacing_y*m_gridSpacing_z;


    const int m_nGridPoints = nElements;
    for(int gId=0; gId<m_nGridPoints; gId++) {
        m_gridPointPorosity[gId] /= nPointsInGridPoint[gId];
        const int np = m_particlesInGridPoint[gId].size();
        const double l_porosity = np*particleVolume/volumegridElement;
        double estimated_nParticles = m_gridPointPorosity[gId]/((double)nElements*solidfraction )*m_nParticles;

        cout << "gId:" << gId << "\tp:" << m_gridPointPorosity[gId] << "\t\tpp:" << l_porosity << "\t\tnp:" << np << "\t\test np:" << estimated_nParticles << endl;
    }

}
//------------------------------------------------------------------------------
void mg::MeshGenerator::computePorosity(mg::Grid &grid, const cube_type & img_data)
{
    const int nElements = grid.nElements();
    vector<double> gridPointPorosity = std::vector<double>(nElements, 0);
    vector<int> nPointsInGridPoint = std::vector<int>(nElements, 0);
    arma::vec3 r = {0,0,0};
    double solidfraction = 0;

    const vector<double> domainLength = grid.getDomainLength();
    const int width = img_data.n_rows;
    const int height = img_data.n_cols;
    const int depth = img_data.n_slices;
    const double dx  = domainLength[0]/width;
    const double dy  = domainLength[1]/height;
    const double dz  = domainLength[2]/depth;

    for(int i=0;i<width;i++) {
        for(int j=0;j<height ;j++) {
            for(int k=0;k<depth;k++) {
                r = {i*dx, j*dy, k*dz};
                const int gId = grid.findGridId(r);
                nPointsInGridPoint[gId] += 1;

                if (img_data(i, j, k) == 0) {
                    gridPointPorosity[gId] += 1;
                    solidfraction  += 1;
                }
            }
        }
    }
    const int nPixels = width*height*depth;

//    solidfraction  /= (nPixels);

    for(int gId=0; gId<nElements; gId++) {
//        gridPointPorosity[gId] /= nPointsInGridPoint[gId];
//        gridPointPorosity[gId] /= solidfraction;
        gridPointPorosity[gId] /= nPixels;
    }
    //-------------------------------
    nPointsInGridPoint = std::vector<int>(nElements, 0);
    vector<double> gridPointDensity = std::vector<double>(nElements, 0); // TMP
    const vector<double> & gridSpacing = grid.getGridspacing();
    const vector<arma::vec3> & gridCenters= grid.getCenters();
    int np = 0;
    solidfraction = 0;

    for(int gId=0; gId<nElements; gId++) {
        const arma::vec3 center = gridCenters[gId];
        double x0 = center[0] - 0.5*gridSpacing[0];
        double x1 = center[0] + 0.5*gridSpacing[0];
        double y0 = center[1] - 0.5*gridSpacing[1];
        double y1 = center[1] + 0.5*gridSpacing[1];
        double z0 = center[2] - 0.5*gridSpacing[2];
        double z1 = center[2] + 0.5*gridSpacing[2];
        for(double x=x0; x<x1; x += dx) {
            for(double y=y0; y<y1; y += dy) {
                for(double z=z0; z<z1; z += dz) {
                    r = {x, y, z};
                    const int gId = grid.findGridId(r);
                    nPointsInGridPoint[gId] += 1;

                    const int i = x/dx;
                    const int j = y/dy;
                    const int k = z/dz;
                    if (img_data(i, j, k) == 0) {
                        gridPointDensity[gId] += 1;
                        solidfraction  += 1;
                    }

                    np++;
                }
            }
        }
    }

    for(int gId=0; gId<nElements; gId++) {
//        gridPointDensity[gId] /= nPointsInGridPoint[gId];
//        gridPointDensity[gId] /= np;
        gridPointDensity[gId] /= solidfraction;
//        cout << gId << " " << gridPointPorosity[gId] << " " << gridPointDensity[gId] << endl;
    }
    //-------------------------------
    grid.registerGridPropery("porosity", gridPointDensity);
//    grid.registerGridPropery("porosity", gridPointPorosity);
//    grid.registerGridPropery("density", gridPointDensity);
}
//------------------------------------------------------------------------------
