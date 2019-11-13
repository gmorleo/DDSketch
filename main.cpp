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
const float DEFAULT_ALPHA = 0.04;

void printDataset(int n_element) {

    ofstream file;
    file.open("normal_mean_2_stdev_3.csv");

    default_random_engine generator;
    normal_distribution<double> distribution(2,3);

    double item;

    for (int i = 0; i < n_element; i++) {
        item = distribution(generator);
        file << item << ", \n";
    }

    file.close();
}

int main() {

    // init the sketch
    DDS_type *dds = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    DDS_type *dds2 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);


    // Test with normal distribution
    int n_element = 1000000000;
    vector<double> stream;
    vector<double> stream2;

    default_random_engine generator;
    normal_distribution<double> distribution(2,3);
    double item;

    // generate pseudorandom number (item) according to the normal distribution
    // insert item into the stream vector and
    // add it to our DDSketch
    for (int i = 0; i < n_element; i++) {
        item = distribution(generator);
        stream.insert(stream.end(), item);
        DDS_Add(dds, item);
    }
/*
    normal_distribution<double> distribution2(1,4);

    for (int i = 0; i < n_element; i++) {
        item = distribution(generator);
        stream2.insert(stream2.end(), item);
        DDS_Add(dds2, item);
    }

    stream.insert(stream.end(), stream2.begin(), stream2.end());
    DDS_merge(dds,dds2);*/

    // sort the vector
    //sort(stream.begin(),stream.end());

    // Determine the quantile
    float q[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};
    int n_q = 10;

    for ( int i = 0; i < n_q; i++ ) {
        int idx = floor(1+q[i]*(stream.size()-1));
        // determine the correct answer
        // i.e., the number stored at the index (idx-1) in the sorted permutation of the vector
        // note that we are not sorting the vector, we are using the quickselect() algorithm
        // which in C++ is available as std::nth_element
        nth_element(stream.begin(), stream.begin() + (idx-1), stream.end());
        double quantile = DDS_GetQuantile(dds, q[i]);
        double error = abs((quantile-stream[idx-1])/stream[idx-1]);
        cout << "q: " << std::setw(4) << q[i] << " estimate: " << std::setw(10) << quantile << " real: " << std::setw(10) << stream[idx-1] << " error: " << std::setw(10) << error << endl;
    }

    // check the number of element counted in bins map
    DDS_SumBins(dds);
    DDS_PrintCSV(dds);
    cout << "max key: " << DDS_GetRank(dds, dds->bins->rbegin()->first) << endl;

    // now check that delete works
    cout << "Sketch size (number of bins) before delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl;
    for (double item : stream) {
        DDS_Delete(dds, item);
        //DDS_CheckAll(dds, item);
    }


    cout << "Sketch size (number of bins) after delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl;

    // deallocate the sketch data structure
    DDS_Destroy(dds);

    printDataset(n_element);


    return 0;
}


