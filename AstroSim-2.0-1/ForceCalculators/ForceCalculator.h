#pragma once
/*
 * INTEGRATORS.H
 */
#ifndef FORCE_CALCULATOR_H_
#define FORCE_CALCULATOR_H_

#include "utils.h"
#include "../Body/body.h"

class ForceCalculator{
    protected:
        GravitationalSystem s;
    public:
        ForceCalculator(GravitationalSystem& s);
        virtual ~ForceCalculator(){
            //cout << "Destroyed Force" << endl;
        }
        virtual valtype getForce(const int i, const int coordType)= 0;
};


#endif