#include "meshgenerator.h"

//------------------------------------------------------------------------------
mg::MeshGenerator::MeshGenerator()
{
}
//------------------------------------------------------------------------------
void mg::save_xyz(arma::mat &x, std::string base = "data_", int i=0)
{
    string fileName = base + to_string(i) + ".xyz";
    ofstream outStream(fileName.c_str());

    outStream << x.n_cols << endl;
    outStream << "#" << endl;

    for (int i=0; i<x.n_cols;i++)
    {
        const arma::vec2& r = x.col(i);
        outStream << i << "\t" << r[0] << "\t" << r[1] << endl;
    }
    outStream.close();
    cout << fileName << endl;
}
//------------------------------------------------------------------------------
void mg::save_image_and_xyz(Parameters parameters, arma::mat & img_data,
                            arma::mat &x, std::string base = "data_", int i=0)
{
    string fileName = base + to_string(i) + ".pgm";

    int n = parameters.nParticles;
    arma::vec areas = arma::zeros(n);

    int h = img_data.n_rows;
    int w = img_data.n_cols;
    double X_0 = 0;
    double X_1 = (double)w/w;
    double Y_0 = 0;
    double Y_1 = (double)h/w;
    double dx = (X_1 - X_0)/w;
    double dy = (Y_1 - Y_0)/h;
    double _dx = (X_1 - X_0);
    double _dy = (Y_1 - Y_0);

    int resolution_x = X_1*parameters.imageResolution;
    int resolution_y = Y_1*parameters.imageResolution;
    arma::mat image(resolution_y, resolution_x);
    image.zeros();

    int pix_hole = 0;

    //--------------------------------------------------------------------------
    // Creating a voronoi image and computing the areas
    //--------------------------------------------------------------------------
    for (int i=0; i<resolution_x;i++)
    {
        for (int j=0; j<resolution_y;j++)
        {
            arma::vec2 r_img;
            r_img[0] = X_1*i/(resolution_x);
            r_img[1] = Y_1*j/(resolution_y);

            double maxLen = numeric_limits<double>::max();
            int indexMax = -1;

            if(img_data(r_img(1)/dy, r_img(0)/dx) > 0){
                image(j, i) = 0;
                pix_hole++;
                continue;
            }

            // Finding the closest Vornoi point
            for (int k=0;k<n; k++)
            {
                arma::vec2 x_k = r_img - x.col(k);

                if(parameters.periodic_x)
                {
                    if(x_k[0] >= 0.5*_dx)
                        x_k[0] -= _dx;
                    if(x_k[0] <= -0.5*_dx)
                        x_k[0] += _dx;
                }

                double dr_rk = x_k[0]*x_k[0] + x_k[1]*x_k[1];

                if(dr_rk < maxLen)
                {
                    maxLen = dr_rk;
                    indexMax = k;
                }
            }
            image(j, i) = indexMax;
            areas[indexMax] += 1.0;
        }
    }

    // Adding the voronoi centers to the image
    for (int k=0;k<n; k++)
        image(x(1, k)*resolution_y/Y_1, x(0, k)*resolution_x) = 1.0;

    // Saving the voronoi image
    image.save(fileName, arma::pgm_binary);

    //--------------------------------------------------------------------------
    // Saving xyz-file with volume
    //--------------------------------------------------------------------------
    fileName = base + to_string(i) + ".xyz";
    ofstream outStream(fileName.c_str());

    outStream << x.n_cols << endl;
    outStream << "# id x y volume" << endl;
    double dxdy = (X_1 - X_0)*(Y_1 - Y_0);
    double height = 1.0;
    double s_volume = 0;
    double total_pix = resolution_x*resolution_y - pix_hole;
    for (int i=0; i<n;i++)
    {
        double volume = dxdy*height* areas[i]/(total_pix);
        const arma::vec2& r = x.col(i);
        outStream << i << "\t" << r[0] << "\t" << r[1] << "\t" << volume << endl;
        s_volume += volume;
    }
    outStream.close();
    cout << fileName << endl;
}
//------------------------------------------------------------------------------
arma::mat mg::createMeshFromImage(Parameters parameters)
{
    CImg<double> image(parameters.imgPath.c_str());

    int h = image.height();
    int w = image.width();
    arma::mat img_data(h, w);

    for(int j=0;j<w;j++)
    {
        for(int i=0;i<h;i++)
        {
            img_data(i,j) = image(j,i,0,0)/(255.0);
        }
    }

    int n = parameters.nParticles;
    int q = parameters.q;
    int threshold = parameters.threshold;

    double a_1 = parameters.alpha_1;
    double a_2 = parameters.alpha_2;
    double b_1 = parameters.beta_1;
    double b_2 = parameters.beta_2;
    double X_0 = 0;
    double X_1 = (double)w/w;
    double Y_0 = 0;
    double Y_1 = (double)h/w;

    arma::mat x = arma::randu(2,n);
    arma::vec js = arma::ones(n);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::uniform_real_distribution<double> distribution_x(X_0, X_1);
    std::uniform_real_distribution<double> distribution_y(Y_0, Y_1);


    double dx = (X_1 - X_0)/w;
    double dy = (Y_1 - Y_0)/h;
    double _dx = (X_1 - X_0);
    double _dy = (Y_1 - Y_0);
    for(int i=0; i<n; i++)
    {
        arma::vec2 r_i = x.col(i);
        do
        {
            r_i[0] = distribution_x(generator);
            r_i[1] = distribution_y(generator);
    }while(img_data(r_i(1)/dy, r_i(0)/dx) > 0);
        x.col(i) = r_i;
    }

    for (int k=0; k<threshold;k++)
    {
        if(parameters.periodic_x)
        {
            for(int k=0;k<n; k++)
            {
                if(x(0,k) < X_0)
                    x(0,k) += _dx;
                if(x(0,k) >= X_1)
                    x(0,k) -= _dx;
            }
        }

        if(k%20 == 0)
        {
            save_image_and_xyz(parameters, img_data, x, parameters.basePath + "/alg1_", k);
        }

        unordered_map <int, vector <arma::vec2>> neighbours;
        for(int r=0; r<q; r++)
        {
            arma::vec2 y_r;
            double maxLen = numeric_limits<double>::max();
            int indexMax = -1;

            do
            {
                y_r[0] = distribution_x(generator);
                y_r[1] = distribution_y(generator);
            }while(img_data(y_r(1)/dy, y_r(0)/dx) > 0);

            arma::vec2 y_tmp = y_r;

            // Finding the closest voronoi center
            for(int k=0;k<n; k++)
            {
                arma::vec2 y_r_copy = y_r;
                arma::vec2 x_k = y_r - x.col(k);

                if(parameters.periodic_x)
                {
                    if(x_k[0] > 0.5*_dx){
                        x_k[0] -= _dx;
                        y_r_copy[0] -= _dx;
                    }else if(x_k[0] < -0.5*_dx){
                        x_k[0] += _dx;
                        y_r_copy[0] += _dx;
                    }
                }

                double dr_rk = x_k[0]*x_k[0] + x_k[1]*x_k[1];

                if(dr_rk < maxLen)
                {
                    maxLen = dr_rk;
                    indexMax = k;
                    y_tmp = y_r_copy;
                }
            }
            neighbours[indexMax].push_back(y_tmp);
        }

        for(int i=0; i<n; i++)
        {
            vector <arma::vec2> neigh = neighbours[i];
            if(neigh.size()>0)
            {
                arma::vec2 x_i =  x.col(i);
                double j = js(i);
                // Finding the average of the random points
                arma::vec2 u_r = arma::zeros(2);
                for (int l=0; l<neigh.size(); l++)
                {
                    const arma::vec2& yr = neigh[l];
                    u_r += yr;
                }
                u_r /= neigh.size();

                x.col(i) = ((a_1*j + b_1)*x_i + (a_2*j + b_2)*u_r)/(j+1);
                js(i) += 1;
            }
        }
    }

    return x;
}
//------------------------------------------------------------------------------
