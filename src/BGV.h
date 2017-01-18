/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BGV.h
 * Author: jnortiz
 *
 * Created on January 17, 2017, 3:19 PM
 */

#ifndef BGV_H
#define BGV_H

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>

using namespace NTL;

class BGV {
public:
    BGV();
    virtual ~BGV();
    
    void HWT(int* &_v, int h);
    void UniformModQ(vec_ZZ& _v, const ZZ& q);

    vec_ZZ GetQ_l() const {
        return _q_l;
    }
    
private:
    void Compute_q_l();
    vec_ZZ _q_l;
};

#endif /* BGV_H */

