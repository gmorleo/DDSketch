#include <iostream>
#include <math.h>
#include <algorithm>
#include <random>
#include <fstream>
#include <iomanip>
#include "ddsketch.h"

/*!
 * \mainpage Project $(project_name) Lorem ipsum dolor
 */
/// \file

using namespace std;

const int DEFAULT_OFFSET = 1073741824; //2^31/2
const int DEFAULT_BIN_LIMIT = 500;
const float DEFAULT_ALPHA = 0.01;

int printDataset(const string& name, int n_element) {

    // open file for output
    ofstream file;
    file.open(name);
    if (file.fail()) {
        cout << "File not open" << endl;;
        return -1;
    }

    // generate with normal distribution
    default_random_engine generator;
    normal_distribution<double> distribution(2,3);

    double item;

    for (int i = 0; i < n_element; i++) {

        item = distribution(generator);
        file << item << ", \n";

    }

    file.close();

    return 0;
}

int insertNormalDistribution(DDS_type *dds, double* stream, int n_element) {

    // Test with normal distribution
    default_random_engine generator;
    normal_distribution<double> distribution(2,3);
    double item;

    // generate pseudorandom number (item) according to the normal distribution
    // insert item into the stream vector and
    // add it to our DDSketch
    for (int i = 0; i < n_element; i++) {

        item = distribution(generator);
        stream[i] = item;
        int return_value = DDS_Add(dds, item);
        //int return_value = DDS_AddRemapped(dds, item);
        if(return_value < 0){
            return -1;
        }
    }

    return 0;
}

int printQuantile(DDS_type *dds, double* stream, int n_element) {

    // Determine the quantile
    float q[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};
    int n_q = 10;

    for ( int i = 0; i < n_q; i++ ) {

        int idx = floor(1+q[i]*(n_element-1));
        // determine the correct answer
        // i.e., the number stored at the index (idx-1) in the sorted permutation of the vector
        // note that we are not sorting the vector, we are using the quickselect() algorithm
        // which in C++ is available as std::nth_element
        nth_element(stream, stream + (idx-1), stream +  n_element);
        double quantile = DDS_GetQuantile(dds, q[i]);
        if (quantile == numeric_limits<double>::quiet_NaN()) {
            cout << "q must be in the interval [0,1]" << endl;
            return  -2;
        }
        double error = abs((quantile-stream[idx-1])/stream[idx-1]);

        cout << "q: " << std::setw(4) << q[i] << " estimate: " << std::setw(10) << quantile << " real: " << std::setw(10) << stream[idx-1] << " error: " << std::setw(10) << error << endl;
    }

    return  0;
}

int deleteElements(DDS_type* dds, double* stream, int n_element) {

    // now check that delete works
    cout << "Sketch size (number of bins) before delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl;
    for (int i = 0; i < n_element; i++) {
        DDS_Delete(dds, stream[i]);
        //DDS_DeleteCollapseNeighbors(dds, stream[i]);
        //DDS_CheckAll(dds, item);
    }

    cout << "Sketch size (number of bins) after delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl;

    DDS_PrintCSV("bins.csv", dds->bins);

    return 0;
}



int main() {

    // init the sketch
    DDS_type *dds = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    // number of element
    int n_element = pow(10,7);

    // init array for all the elements
    double* stream = NULL;
    stream = new (nothrow) double[n_element];
    if(!stream){
        cout << "Not enough memory" << endl;;
        return -2;
    }

    insertNormalDistribution(dds, stream, n_element);
    printQuantile(dds, stream, n_element);
    deleteElements(dds, stream, n_element);

    // print dataset on a file
    //printDataset("normal.csv", n_element);

    // deallocate the sketch data structure
    DDS_Destroy(dds);

    return 0;
}


