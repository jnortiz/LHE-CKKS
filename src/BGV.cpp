/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BGV.cpp
 * Author: jnortiz
 * 
 * Created on January 17, 2017, 3:19 PM
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "BGV.h"
#include "params.h"
#include <unistd.h>
#include <sys/syscall.h>
#include <linux/random.h>
#include <cinttypes>
#include <cmath>

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_RR.h>

using namespace std;

BGV::BGV() {
    
    this->_q_l.SetLength(L_def+1);
    this->Compute_q_l();
    
}

BGV::~BGV() {
}

void BGV::Compute_q_l(void) {
    
    /* Invalid position */
    this->_q_l[0] = 0;
    
    /* Following definition on page 7 of [CKKS, 2016] */
    for(int i = 1; i <= L_def; i++) {    
        this->_q_l[i] = NTL::to_ZZ(pow((double) p_def, (double) i)*q0_def);
    }

}
/**
 * 
 * @param h - exactly the numbers of nonzero entries in the vector _v
 * @return _v - a vector with exactly h nonzero entries
 */
void BGV::HWT(int* &_v, int h) {
    /* It's desirable that h is approximately N/2 considering that
     the distribution of 0s and 1s is similar in a uniform distribution */
    
    _v = new int[n_def];
    
    int hw, it_v;    
    size_t l = n_def/64;
    uint64_t *_buf = new uint64_t[l];    
    uint64_t *_signal = new uint64_t[l];
    
    do {
                
        syscall(SYS_getrandom, _buf, l*8, (unsigned int) 1);
        syscall(SYS_getrandom, _signal, l*8, (unsigned int) 1);
        
        it_v = 0;
        hw = 0;
        for(int i = 0; i < l; i++) { /* For each 64-bit word */
            for(int j = 0; j < 64; j++) {
                /* Expand the word as a binary vector */
                _v[it_v] = _buf[i] & 0x1;
                /* Compute the Hamming weight of the entire N-bit vector */
                hw += _v[it_v];
                
                /* Determine if a nonzero entry is positive or negative by using 
                 a random word of the _signal vector */
                if(_signal[i] & 0x1) {
                    _v[it_v] *= -1;
                }
                
                _signal[i] >>= 1;
                _buf[i] >>= 1;
                it_v++;
            }
        }
        
    } while(hw != h);
   
    delete [] _buf;
    delete [] _signal;
    _buf = NULL;
    _signal = NULL;

}

//
/*
 * Uniform distribution over Z_q^N 
 * @param _v
 * @param q 
 */
void BGV::UniformModQ(vec_ZZ& _v, const ZZ& q) {
    
    _v.SetLength(n_def);
    ZZ rnd;
    
    for(int i = 0; i < n_def; i++) {
        rnd = RandomBnd(q);
        /* Rescaling into the interval (-q/2, q/2] */
        _v[i] = (rnd - to_ZZ(q/2));
    }
    
}

/* If the bit is zero, the output becomes "a" */
int Select(int a, int b, unsigned bit) {
    
    unsigned mask;
    int output;
    
    mask = -bit;
    output = mask & (a ^ b);
    output = output ^ a;
    
    return output;
    
}//end-Select()

RR Probability(RR x, RR sigma, RR c) {
    
    RR S = sigma*sqrt(2*ComputePi_RR());
    RR overS = 1/S;
    
    if(x == to_RR(0))
        return overS;
    
    return overS*exp(-(power((x-c)/sigma, 2))/2.0);
    
}//end-Probability()

int KnuthYao(const Vec< Vec<int> >& P, const Vec<int>& begin, int tailcut, RR sigma, RR c, int q) {

    int bound, center, col, d, invalidSample, pNumRows, pNumCols, S, signal;
    unsigned enable, hit;
    unsigned long r;
    
    bound = tailcut*to_int(sigma);
    center = to_int(c);
    d = 0;
    hit = 0;
    signal = 1 - 2*RandomBits_long(1);
    invalidSample = bound+1;
    pNumRows = P.length();
    pNumCols = P[0].length();    
    
    Vec<int> _randomBits;
    _randomBits.SetLength(pNumRows);
    
    int i, index, j, length;
    length = sizeof(unsigned long)*8; 
    
    index = 0;
    for(i = 0; i < (pNumRows/length+1); i++) {
        r = RandomWord();
        for(j = 0; j < length && index < pNumRows; j++, r >>= 1)
            _randomBits[index++] = (r & 1);
    }//end-for

    S = 0;    
    for(int row = 0; row < pNumRows; row++) {
        
        d = 2*d + _randomBits[row]; // Distance calculus
        
        for(col = begin[row]; col < pNumCols; col++) {
            
            d = d - P[row][col];
            
            enable = (unsigned)(d + 1); // "enable" turns 0 iff d = -1
            enable = (1 ^ ((enable | -enable) >> 31)) & 1; // "enable" turns 1 iff "enable" was 0
             
            /* When enable&!hit becomes 1, "col" is added to "S";
             * e.g. enable = 1 and hit = 0 */
            S += Select(invalidSample, col, (enable & !hit));
            hit += (enable & !hit);
                            
        }//end-for
        
    }//end-for
    
    /* Note: the "col" value is in [0, bound]. So, the invalid sample must be 
     * greater than bound. */
    S %= invalidSample;
    S = S - bound + center;
    S *= signal;
    
    return (int)(S - q/2); /* Z_q = Z intersection with (-q/2, q/2] */
    
}//end-Knuth-Yao()

/* This method build the probability matrix for samples in the range 
 * [-tailcut*\floor(sigma), +tailcut*\floor(sigma)] */
void BuildProbabilityMatrix(Vec< Vec<int> >& P, Vec<int>& begin, int precision, int tailcut, RR sigma, RR c) {
    
    RR::SetPrecision(to_long(precision));

    Vec< Vec<int> > __auxP;
    Vec<int> _auxBegin;
    
    // The random variable consists of elements in [c-tailcut*sigma, c+tailcut*sigma]
    int i, j, bound, pNumCols, pNumRows, x;
    vec_RR probOfX;
    RR pow;
    
    bound = to_int(to_RR(tailcut)*sigma);
    
    probOfX.SetLength(bound+1);
       
    __auxP.SetLength(precision);
    for(i = 0; i < __auxP.length(); i++)
        __auxP[i].SetLength(bound+1);

    for(x = bound; x > 0; x--)
        probOfX[bound-x] = Probability(to_RR(x) + c, sigma, c);
    div(probOfX[bound], Probability(to_RR(0) + c, sigma, c), to_RR(2));
    
    i = -1;
    for(j = 0; j < precision; j++) {
        pow = power2_RR(i--); // 2^{i}
        for(x = bound; x >= 0; x--) {
            __auxP[j][bound-x] = 0;                
            if(probOfX[bound-x] >= pow) {
                __auxP[j][bound-x] = 1;
                probOfX[bound-x] -= pow;
            }//end-if
        }//end-for
    }//end-while
    
    P = __auxP;
    
    pNumCols = P[0].length();
    pNumRows = P.length();
    
    _auxBegin.SetLength(pNumRows);
    
    // Computing in which position the non-zero values in P start and end 
    for(i = 0; i < pNumRows; i++) {
        
        _auxBegin[i] = pNumCols-1;
        
        for(j = 0; j < pNumCols; j++)
            if(P[i][j] == 1) {
                _auxBegin[i] = j;
                break;
            }//end-if
        
    }//end-for
    
    begin = _auxBegin;
                
}//end-BuildProbabilityMatrix()
