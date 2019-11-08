#include <iostream>
#include <math.h>
#include <algorithm>
#include <random>
#include "ddsketch.h"

/*!
 * \mainpage Project $(project_name) Lorem ipsum dolor
 */
/// \file

using namespace std;

const int DEFAULT_OFFSET = 1073741824; //2^31/2
const int DEFAULT_BIN_LIMIT = 500;
const float DEFAULT_ALPHA = 0.01;


int main() {

    // init the sketch
    DDS_type *dds = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);


    // Test with normal distribution
    int n_element = 10000000;
    vector<double> stream;

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

    // sort the vector
    //sort(stream.begin(),stream.end());

    // determine the 0.7 quantile
    float q = 0.7;
    int idx = floor(1+q*(stream.size()-1));

    // determine the correct answer
    // i.e., the number stored at the index (idx-1) in the sorted permutation of the vector
    // note that we are not sorting the vector, we are using the quickselect() algorithm
    // which in C++ is available as std::nth_element
    nth_element(stream.begin(), stream.begin() + (idx-1), stream.end());

    double quantile = DDS_GetQuantile(dds, q);
    double error = abs((quantile-stream[idx-1])/stream[idx-1]);

    cout << "Result for q = " << q << endl;
    cout << "Real: " << stream[idx-1] << " Estimation: " << quantile << " Error: " << error << endl;


    // check the number of element counted in bins map
    DDS_SumBins(dds);
    DDS_PrintCSV(dds);

    // now check that delete works
    cout << "Sketch size (number of bins) before delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl;
    for (double item : stream) {
        //DDS_Delete(dds, item);
        DDS_CheckAll(dds, item);
    }


    cout << "Sketch size (number of bins) after delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl;

    // deallocate the sketch data structure
    DDS_Destroy(dds);


    return 0;
}


