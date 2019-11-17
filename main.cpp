/********************************************************************
 DDSketch

 An algorithm for tracking quantiles in data streams

 Charles Masson, Jee E. Rim, and Homin K. Lee. 2019. DDSketch: a fast and fully-mergeable quantile sketch with relative-error guarantees. Proc. VLDB Endow. 12, 12 (August 2019), 2195-2205. DOI: https://doi.org/10.14778/3352063.3352135

 This implementation by
 by Giuseppe Morleo
 University of Salento, Italy

 *********************************************************************/
/*!
 * \mainpage Project $(project_name) DDSketch
 */
/// \file

#include <iostream>
#include <math.h>
#include <algorithm>
#include <random>
#include <fstream>
#include <iomanip>
#include "ddsketch.h"



using namespace std;

const int DEFAULT_OFFSET = 1073741824; //2^31/2
const int DEFAULT_BIN_LIMIT = 500;
const double DEFAULT_ALPHA = 0.008;

/**
 * @brief               This function generates a csv file of (n_element) based on the normal distribution
 * @param name          Output file
 * @param n_element     Number of element in the dataset
 * @return              0: success; \n-1: file opening failed;
 */
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

/**
 * @brief               This fucntion inserts (n_element) to the sketch according the normal distribution
 * @param dds           Parameters of the sketch
 * @param stream        Vector that contains all the real values inserted
 * @param n_element     Number of element
 * @return              0: success; \n-1: error;
 */
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

/**
 * @brief               This function computes the quantile
 * @param dds           Parameters of the sketch
 * @param stream        Vector that contains all the real values inserted
 * @param n_element     Number of element
 * @return              0: success; \n-2: error;
 */
int printQuantile(DDS_type *dds, double* stream, int n_element) {

    // Determine the quantile
    float q[] = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};
    int n_q = 11;

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

/**
 * @brief               This function checks if all elements in the stream have a corresponding bucket in the sketch
 * @param dds           Parameters of the sketch
 * @param stream        Vector that contains all the real values inserted
 * @param n_element     Number of element
 * @return              0: success
 */
int deleteElements(DDS_type* dds, double* stream, int n_element) {

    // now check that delete works
    cout << "Sketch size (number of bins) before delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl;

    for (int i = 0; i < n_element; i++) {
        //DDS_DeleteCollapseNeighbors(dds, stream[i]);
        //DDS_CheckAll(dds, item);
        int return_value = DDS_Delete(dds, stream[i]);
        if ( return_value<0 ) {
            cout << "Key associated to the value " << stream[i] << " not found!" << endl;
        }
    }

    cout << "Sketch size (number of bins) after delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl;

    return 0;
}


/**
 *
 * @return
 */
int main() {

    // init the sketch
    DDS_type *dds = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    // number of element
    int n_element = pow(10,8);

    // init array for all the elements
    double* stream = nullptr;
    stream = new (nothrow) double[n_element];
    if(!stream){
        cout << "Not enough memory" << endl;
        DDS_Destroy(dds);
        return -2;
    }

    insertNormalDistribution(dds, stream, n_element);
    printQuantile(dds, stream, n_element);
    deleteElements(dds, stream, n_element);

    // print dataset on a file
    //printDataset("normal.csv", n_element);

    // deallocate the sketch data structure
    DDS_Destroy(dds);

    delete[] stream, stream = nullptr;

    return 0;
}


