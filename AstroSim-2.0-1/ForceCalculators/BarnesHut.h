#pragma once
/*
 * INTEGRATORS.H
 */
#ifndef BARNES_HUT_H_
#define BARNES_HUT_H_

#include "ForceCalculator.h"

class BarnesHut : public ForceCalculator{
    private:
        vector<vector<long double>> forces;
    public:
        BarnesHut(GravitationalSystem& s);
        ~BarnesHut(){
            //cout << "Destroyed BanresHut" << endl;
        }
        valtype getForce (const int i, const int coordType) override;
};


#endif