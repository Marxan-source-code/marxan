#pragma once

namespace marxan {
    using namespace std;

    typedef struct spu
    {
        double amount; // amount of species in pu
        double prob; // optional field that does ???
        int spindex; // index/id of species
    } spu;

    //Separate output fields for multithreading 
    typedef struct spu_out
    {
        int clump;
        spu_out() { clump = 0; }
    } spu_out;

    typedef struct spusporder
    {
        double amount;
        int puindex;
    } spusporder;

    // type definitions for original Ian Ball data structures

    typedef struct spustuff
    {
        int id;
        int status;
        double xloc, yloc;
        double cost;
        double prob;
        int richness, offset, probrichness, proboffset;
    } spustuff;

} // namespace marxan