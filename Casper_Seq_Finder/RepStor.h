//
//  RepStor.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 10/23/20.
//  Copyright © 2020 Brian Mendoza. All rights reserved.
//

#ifndef RepStor_h
#define RepStor_h

#include <iostream>
#include <vector>
#include <set>
#include "gRNA.h"


class RepStor {
public:
    RepStor();
    ~RepStor();
    
public:
    void addNewRepeat(gRNA*);  //adds a new gRNA object to the storage container
    
private:
    gRNA* dual [2];
    gRNA* tohundred [98];
    
};

#endif /* RepStor_h */
