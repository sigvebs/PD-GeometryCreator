#include "grid.h"

#include <math.h>

namespace mg
{
//------------------------------------------------------------------------------
Grid::Grid(double x0, double x1, double y0, double y1, double z0, double z1, double m_approxSpacing):
    m_domain({x0, x1, y0, y1, z0, z1}),
    m_domainLength({x1 - x0, y1 - y0, z1 - z0}),
    m_approxSpacing(m_approxSpacing)
{

}
//------------------------------------------------------------------------------
void Grid::createGrid()
{
    m_nGrid = {0,0,0};
    m_gridSpacing = {0,0,0};
    m_nGridElements = 1;
    for(int d=0; d<M_DIM; d++) {
        m_nGrid[d] = floor(m_domainLength[d]/m_approxSpacing) > 0 ? floor(m_domainLength[d]/m_approxSpacing) : 1;
        m_gridSpacing[d] = m_domainLength[d]/m_nGrid[d];
        m_nGridElements *= m_nGrid[d];
    }

    m_centers = vector<arma::vec3>(m_nGridElements, arma::vec3({0,0,0}));
    // Setting the grid centers
    for(int i=0; i<m_nGrid[0]; i++) {
        for(int j=0; j<m_nGrid[1]; j++) {
            for(int k=0; k<m_nGrid[2]; k++) {
                const int id_ijk = k + m_nGrid[2]*j + m_nGrid[1]*m_nGrid[2]*i;
                arma::vec3 center = {(i+0.5)*m_gridSpacing[0], (j+0.5)*m_gridSpacing[1], (k+0.5)*m_gridSpacing[2]};
                m_centers[id_ijk] = center;
            }
        }
    }
}
//------------------------------------------------------------------------------
int Grid::findGridId(const arma::vec3 &r)
{
    std::vector<int> id_xyz = {r(0)/m_gridSpacing[0], r(1)/m_gridSpacing[1], r(2)/m_gridSpacing[2]};

    // Boundary checks
    for(int d=0; d<M_DIM; d++) {

        if(id_xyz[d] > m_nGrid[d])
            id_xyz[d] = m_nGrid[d] - 1;
        else if(id_xyz[d] < 0)
            id_xyz[d] = 0;
    }

    return id_xyz[2] + m_nGrid[2]*id_xyz[1] + m_nGrid[1]*m_nGrid[2]*id_xyz[0];
}
//------------------------------------------------------------------------------
void Grid::registerGridPropery(const string &str)
{
    m_gridProperties[str] = vector<double>(m_nGridElements, 0);
}
//------------------------------------------------------------------------------
void Grid::registerGridPropery(const string &str, vector<double> data)
{
    m_gridProperties[str] = data;
}
//------------------------------------------------------------------------------
const vector<double> &Grid::getGridProperties(const string &str)
{
    if(m_gridProperties.count(str) == 0) {
        cerr << "WARNING: Parameter '" << str
             << "' is not registered." << endl;
    }

    return m_gridProperties[str];
}
//------------------------------------------------------------------------------
int Grid::nElements()
{
    return m_nGridElements;
}
//------------------------------------------------------------------------------
const vector<double> &Grid::getDomainLength()
{
    return m_domainLength;
}
//------------------------------------------------------------------------------
const vector<arma::vec3> &Grid::getCenters()
{
    return m_centers;
}
//------------------------------------------------------------------------------
const vector<double> &Grid::getGridspacing()
{
    return m_gridSpacing;
}
//------------------------------------------------------------------------------
}
