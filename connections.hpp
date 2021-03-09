#pragma once
#include <vector>

namespace marxan {
    using namespace std;

    typedef struct scost
    {
        double total;
        int pus;
        double connection;
        int missing;
        double penalty;
        double cost;
        double threshpen;
        double shortfall;
        double probability1D;
        double probability2D;
    } scost;

    /* Connectivity Structure. Fixed connectivity number.*/
    typedef struct sneighbour
    {
        int nbr; // puid of neighbour
        double cost;
        int connectionorigon; // for asymmetric connections, whether connection starts from this node.
        sneighbour(int nbr, double cost, int connectionorigon) { this->nbr = nbr; this->cost = cost; this->connectionorigon = connectionorigon; }
    } sneighbour;

    typedef struct sconnections
    {
        vector<sneighbour> first;
        double fixedcost;
        int nbrno;
    } sconnections;

} // namespace marxan