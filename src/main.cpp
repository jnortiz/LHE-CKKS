/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on January 17, 2017, 2:58 PM
 */

#include <cstdlib>
#include <iostream>
#include "params.h"

#include "BGV.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    BGV *inst = new BGV();
    int *_v;
    
    inst->HWT(_v, 32);
    
    for(int i = 0; i < N; i++) {
        cout << _v[i] << " ";
    }
    
    return 0;
}

