//
// Created by giuseppe on 25/11/19.
//

#include <iostream>

const int GENERIC_ERROR = -1;
const int MEMORY_ERROR = -2;
const int FILE_ERROR = -3;
const int SKETCH_ERROR = -4;
const int MERGE_ERROR = -5;
const int QUANTILE_ERROR = -6;

using namespace std;

extern int printError(int error, string nameFunction);
