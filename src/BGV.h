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

class BGV {
public:
    BGV();
    BGV(const BGV& orig);
    virtual ~BGV();
    void HWT(int* &_v, int h); 
private:
};

#endif /* BGV_H */

