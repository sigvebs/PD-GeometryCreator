#ifndef GRID_H
#define GRID_H

#include<vector>
#include<armadillo>
#include <unordered_map>

using namespace std;

namespace mg
{
static const int M_DIM = 3;
//------------------------------------------------------------------------------
class Grid
{
protected:
    vector<double> m_domain;
    vector<double> m_domainLength;
    double m_approxSpacing;
    vector<double> m_gridSpacing;
    vector<int> m_nGrid; // n points in direction
    unordered_map<string, vector<double>> m_gridProperties;
    int m_nGridElements;
    vector<arma::vec3> m_centers;

public:
    Grid(double x0, double x1, double y0, double y1, double z0, double z1, double m_approxSpacing);
    void createGrid();
    int findGridId(const arma::vec3 & r_i);
    void registerGridPropery(const string & str);
    void registerGridPropery(const string & str, vector<double> data);
    const vector<double> & getGridProperties(const string & str);
    int nElements();
    const vector<double> &getDomainLength();
    const vector<arma::vec3> &getCenters();
    const vector<double> &getGridspacing();
};
//------------------------------------------------------------------------------
}
#endif // GRID_H
