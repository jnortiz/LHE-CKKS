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

using namespace std;

BGV::BGV() {
}

BGV::BGV(const BGV& orig) {
}

BGV::~BGV() {
}

/**
 * 
 * @param h - exactly the numbers of nonzero entries in the vector _v
 * @return _v
 */
void BGV::HWT(int* &_v, int h) {

    _v = new int[N];
    
    int hw, it_v;    
    size_t l = N/64;
    uint64_t *buf = new uint64_t[l];    
    int64_t *signal = new int64_t[l];
    
    syscall(SYS_getrandom, signal, l*8, (unsigned int) 1);
    for(int i = 0; i < l; i++) { /* For each 64-bit word */
        cout << signal[i] << endl;
    }
    
    do {
                
        syscall(SYS_getrandom, buf, l*8, (unsigned int) 1);
        syscall(SYS_getrandom, signal, l*8, (unsigned int) 1);
        
        it_v = 0;
        hw = 0;
        for(int i = 0; i < l; i++) { /* For each 64-bit word */
            for(int j = 0; j < 64; j++) {
                /* Expand the word as a binary vector */
                _v[it_v] = buf[i] & 0x1;
                /* Compute the Hamming weight of the entire vector */
                hw += _v[it_v];
                
                if(signal[i] & 0x1) {
                    cout << _v[it_v] << " ";
                    _v[it_v] *= -1;
                }
                
//                cout << _v[it_v] << " ";
                signal[i] >>= 1;
                buf[i] >>= 1;
                it_v++;
            }
        }
        cout << endl << "hw: " << hw << endl;
        
    } while(hw != h);
   
    delete [] buf;
    buf = NULL;

}

