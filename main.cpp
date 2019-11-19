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
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define RESET   "\033[0m"

using namespace std;

const double NORMAL_MEAN = 2;
const double NORMAL_STDDEV = 3;
const double EXPONENTIAL_LAMBDA = 17;
const double REAL_UNIFORM_A = -50;
const double REAL_UNIFORM_B = 50;
const double GAMMA_ALPHA = 1;
const double GAMMA_BETA = 3;
const double CHI_SQUARED_N = 10;
const double WEIBULL_A = -4;
const double WEIBULL_B = 4;

const int DEFAULT_OFFSET = 1073741824; //2^31/2
const int DEFAULT_BIN_LIMIT = 500;
const double DEFAULT_ALPHA = 0.008;

/**
 * \brief                   This function computes the dimension of the dataset
 * @param name_file         Name of dataset
 * @return                  Return the number of element in the dataset(row)
 */
int getDatasetSize(const string &name_file);

/**
 * \brief                   This function loads the dataset into an array
 * @param name_file         Name of dataset
 * @param dataset           Array
 * @return                  An array containing the whole dataset
 */
int loadDataset(const string &name_file, double *dataset, int start);

/**
 * @brief                   This function returns the distribution name
 * @param distribution      Code of distribution
 * @return                  Distribution name
 */
string getDistributionName(int distribution);

/**
 * @brief                   This function generates a csv file of (n_element) based on the normal distribution
 * @param name              Output file
 * @param n_element         Number of element in the dataset
 * @return                  0 success; \n-1 file opening failed;
 */
int printDataset(const string& name, int n_element, int distribution);

/**
 * @brief                   This function inserts (n_element) to the sketch according to the selected distribution and collapse type
 * @param dds               Parameters of the sketch
 * @param stream            Vector that contains all the real values inserted
 * @param n_element         Number of element
 * @return                  0 success; \n-1 error;
 */

/**
 * @brief                   This function inserts (n_element) to the sketch according to the selected distribution and collapse type
 * @param dds               The sketch
 * @param stream            Vector that contains all the real values inserted
 * @param stream_start      Initial index from which insert the elements in the vector
 * @param n_element         Number of element
 * @param distribution      Distribution type: 1 normal, 2 exponential, 3 real uniform, 4 gamma, 5 chi squared, 6 weibull
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
int insertRandom(DDS_type *dds, double* stream, int stream_start, int n_element, int distribution, int collapseType);

/**
 * @brief                   This function inserts elements from index (start) to index (stop) from dataset to the sketch according to the selected collapse type
 * @param dds               The sketch
 * @param start             Dataset initial index
 * @param stop              Dataset final index
 * @param dataset           Dataset
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
int insertFromDataset(DDS_type* dds, int start, int stop, double* dataset, int collapseType);

/**
 * @brief                   This function merge two sketches based on the used collapse type
 * @param dds1              First sketch
 * @param dds2              Second sketch
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucke
 * @return                  0 success -1 error
 */
int merge(DDS_type* dds1, DDS_type* dds2, int collapseType);

/**
 * @brief                   This function tests the sketch, entering (n_element) based on the selected distribution and using the selected type of collapsing
 * @param n_element         Number of element
 * @param distribution      Distribution type: 1 normal, 2 exponential, 3 real uniform, 4 gamma, 5 chi squared, 6 weibull
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
int testWithRandomValue(int n_element, int distribution, int collapseType);

/**
 * @brief                   This function tests the merge of two sketches, entering (n_element1) in the first sketch and
 *                          (n_element2) in the second sketch based on the selected distributions and using the selected type of collapsing
 * @param n_element1        Number of element first distribution
 * @param n_element2        Number of element second distribution
 * @param distribution1     First distribution type: 1 normal, 2 exponential, 3 real uniform, 4 gamma, 5 chi squared, 6 weibull
 * @param distribution2     Second distribution type: 1 normal, 2 exponential, 3 real uniform, 4 gamma, 5 chi squared, 6 weibull
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
int testMergeWithRandomValue(int n_element1, int n_element2, int distribution1, int distribution2, int collapseType);

/**
 * @brief                   This function tests the merge of two sketches, entering (dataset size/2) in the first sketch
 *                          and (dataset size/2) in the second sketch from a dataset and using the selected type of collapsing
 * @param dataset_name      Dataset name
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
int testMergeFromDataset(const string& dataset_name, int collapseType);

/**
 * @brief                   This function tests the merge of two sketches, entering (n_element1) in the first sketch and
 *                          (n_element2) in the second sketch from the dataset and using the selected type of collapsing
 * @param dataset1_name     Name first dataset
 * @param dataset2_name     Name second dataset
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
int testMergeFromTwoDataset(const string& dataset1_name, const string& dataset2_name, int collapseType);

/**
 * @brief                   This function checks if all elements in the stream have a corresponding bucket in the sketch
 * @param dds               Parameters of the sketch
 * @param stream            Vector that contains all the real values inserted
 * @param n_element         Number of element
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
int deleteElements(DDS_type* dds, double* stream, int n_element, int collapseType);

/**
 * @brief               This function computes the quantile
 * @param dds           Parameters of the sketch
 * @param stream        Vector that contains all the real values inserted
 * @param n_element     Number of element
 * @return              0: success; \n-2: error;
 */
int testQuantile(DDS_type *dds, double* stream, int n_element);

/**
 *
 * @return
 */
int main() {

    /// number of element
    int n_element = pow(10,8);

    /// Test with random value
    /// Distribution: 1 Normal, 2 Exponential, 3 Real Uniform,4 Gamma, 5 Chi Squared, 6 Weibull
    /// Collapse type: 1 Collapse with gamma^2, 2 Collapse with last buckets, 3 Collapse with first buckets

    testWithRandomValue(n_element, 1, 1);
    //testMergeWithRandomValue(n_element,n_element,1,1,1);

    /// Test merge function
    /// Collapse type: 1 Collapse with gamma^2, 2 Collapse with last buckets, 3 Collapse with first bucket

    //testMergeFromDataset("normal.csv", 1);
    //testMergeFromTwoDataset("normal.csv", "exponential.csv", 1);

    /// Print dataset on a file
    /// Distribution: 1 Normal, 2 Exponential, 3 Real Uniform,4 Gamma, 5 Chi Squared, 6 Weibull

    //printDataset("normal.csv", 10000000, 1);

    return 0;
}

string getDistributionName(int distribution) {
    switch (distribution) {
        case 1:
            return "normal";
        case 2:
            return "exponential";
        case 3:
            return "uniform real";
        case 4:
            return "gamma";
        case 5:
            return "chi sqared";
        case 6:
            return "weibull";
        default:
            return "Wrong distribution code";
    }
}

int printDataset(const string& name, int n_element, int distribution) {

    // open file for output
    ofstream file;
    file.open(name);
    if (file.fail()) {
        cout << "File not open" << endl;;
        return -1;
    }

    // Init the distribution
    default_random_engine generator;
    normal_distribution<double> normal(NORMAL_MEAN,NORMAL_STDDEV);
    exponential_distribution<double> exponential(EXPONENTIAL_LAMBDA);
    uniform_real_distribution<double> uniform_real(REAL_UNIFORM_A,REAL_UNIFORM_B);
    gamma_distribution<double> gamma(GAMMA_ALPHA,GAMMA_BETA);
    chi_squared_distribution<double> chi_squared(CHI_SQUARED_N);
    weibull_distribution<double> weibull(WEIBULL_A,WEIBULL_B);

    double item;

    for (int i = 0; i < n_element; i++) {

        // generate pseudorandom number (item) according to the select distribution
        switch (distribution) {
            case 1:
                item = normal(generator);
                break;
            case 2:
                item = exponential(generator);
                break;
            case 3:
                item = uniform_real(generator);
                break;
            case 4:
                item = gamma(generator);
                break;
            case 5:
                item = chi_squared(generator);
                break;
            case 6:
                item = weibull(generator);
                break;
            default:
                return -1;
        }

        // insert number (item) to the file
        file << item << ", \n";

    }

    file.close();

    return 0;
}

int insertRandom(DDS_type *dds, double* stream, int stream_start, int n_element, int distribution, int collapseType) {

    // Init the distribution
    default_random_engine generator;
    normal_distribution<double> normal(NORMAL_MEAN,NORMAL_STDDEV);
    exponential_distribution<double> exponential(EXPONENTIAL_LAMBDA);
    uniform_real_distribution<double> uniform_real(REAL_UNIFORM_A,REAL_UNIFORM_B);
    gamma_distribution<double> gamma(GAMMA_ALPHA,GAMMA_BETA);
    chi_squared_distribution<double> chi_squared(CHI_SQUARED_N);
    weibull_distribution<double> weibull(WEIBULL_A,WEIBULL_B);

    double item;

    // generate pseudorandom number (item) according to the select distribution
    // insert item into the stream vector and
    // add it to our DDSketch
    for (int i = 0; i < n_element; i++) {

        switch (distribution) {
            case 1:
                item = normal(generator);
                break;
            case 2:
                item = exponential(generator);
                break;
            case 3:
                item = uniform_real(generator);
                break;
            case 4:
                item = gamma(generator);
                break;
            case 5:
                item = chi_squared(generator);
                break;
            case 6:
                item = weibull(generator);
                break;
            default:
                cout << "Wrong distribution number" << endl;
                return -1;
        }

        stream[i+stream_start] = item;

        switch (collapseType) {
            case 1:
                DDS_AddCollapse(dds, item);
                break;
            case 2:
                DDS_AddCollapseLastBucket(dds, item);
                break;
            case 3:
                DDS_AddCollapseFirstBucket(dds, item);
                break;
            default:
                cout << "Wrong collapse type number" << endl;
                return -1;
        }

    }

    return 0;
}

int insertFromDataset(DDS_type* dds, int start, int stop, double* dataset, int collapseType) {

    // insert item into  our DDSketch
    if ( collapseType == 1 ) {
        for (int i = start; i < stop; i++) {
            DDS_AddCollapse(dds, dataset[i]);
        }
    } else if ( collapseType == 2 ) {
        for (int i = start; i < stop; i++) {
            DDS_AddCollapseLastBucket(dds, dataset[i]);
        }
    } else if ( collapseType == 3 ) {
        for (int i = start; i < stop; i++) {
            DDS_AddCollapseFirstBucket(dds, dataset[i]);
        }
    } else {
        cout << "Wrong collapse type number" << endl;
        return -1;
    }

    return 0;
}

int merge(DDS_type* dds1, DDS_type* dds2, int collapseType) {

    // selects the merge function based on the selected collapse method
    switch (collapseType) {
        case 1:
            DDS_MergeCollapse(dds1, dds2);
            break;
        case 2:
            DDS_MergeCollapseLastBucket(dds1, dds2);
            break;
        case 3:
            DDS_MergeCollapseFirstBucket(dds1, dds2);
            break;
        default:
            return -1;
    }

    return 0;
}

int getDatasetSize(const string &name_file) {

    ifstream inputFile(name_file);
    if (inputFile.fail()) {
        cout << "Error " << name_file << " not opened" << endl;
        return -5;
    }
    string line;

    long rows = 0;

    while (getline(inputFile, line))
        rows++;

    return rows;
}

int loadDataset(const string &name_file, double *dataset, int start) {

    ifstream inputFile(name_file);
    if (inputFile.fail()) {
        cout << "Error " << name_file << " not opened" << endl;
        return -5;
    }

    string line;

    int row = start;

    while (getline(inputFile, line)) {
        dataset[row] = stod(line);
        row++;
    }

    return 0;
}

int testWithRandomValue(int n_element, int distribution, int collapseType) {

    cout << endl << BOLDRED << "Test with distribution: " << getDistributionName(distribution) << " initial alpha = " << DEFAULT_ALPHA << " bin limit = " << DEFAULT_BIN_LIMIT << RESET << endl;

    // init the sketch
    DDS_type* dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    // init array for all the elements
    double* stream = nullptr;
    stream = new (nothrow) double[n_element];
    if(!stream){
        cout << "Not enough memory" << endl;
        DDS_Destroy(dds1);
        return -2;
    }

    // Test with random value
    insertRandom(dds1, stream, 0, n_element, distribution, collapseType);
    testQuantile(dds1, stream, n_element);
    deleteElements(dds1, stream, n_element, collapseType);

    DDS_Destroy(dds1);
    delete[] stream, stream = nullptr;

    return 0;
}

int testMergeWithRandomValue(int n_element1, int n_element2, int distribution1, int distribution2, int collapseType) {


    cout << endl << BOLDRED << "Test with distribution: " << getDistributionName(distribution1) << " and " << getDistributionName(distribution2)  << " initial alpha = " << DEFAULT_ALPHA << " bin limit = " << DEFAULT_BIN_LIMIT << RESET << endl;

    // init the sketch
    DDS_type* dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    DDS_type* dds2 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    // init array for all the elements
    double* stream = nullptr;
    stream = new (nothrow) double[n_element1+n_element2];
    if(!stream){
        cout << "Not enough memory" << endl;
        DDS_Destroy(dds1);
        return -2;
    }

    // Test merge with two sketches
    cout << endl << BOLDRED << "Test with two sketches" << RESET << endl;

    cout << "first sketch: " << endl;
    insertRandom(dds1, stream, 0, n_element1, distribution1, collapseType);
    cout << "second sketch: " << endl;
    insertRandom(dds2, stream, n_element1, n_element2, distribution2, collapseType);
    merge(dds1, dds2, collapseType);
    testQuantile(dds1, stream, n_element1 + n_element2);
    deleteElements(dds1, stream, n_element1 + n_element2, collapseType);

    // reset dds1
    DDS_Destroy(dds1);
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    cout << BOLDRED << "Test with one sketch" << RESET << endl;

    // Test result with only one sketch
    insertFromDataset(dds1, 0, n_element1 + n_element2, stream, collapseType);
    testQuantile(dds1, stream, n_element1 + n_element2);
    deleteElements(dds1, stream, n_element1 + n_element2, collapseType);

    // deallocate the sketch data structure
    DDS_Destroy(dds1);
    DDS_Destroy(dds2);

    delete[] stream, stream = nullptr;
    return 0;
}

int testMergeFromDataset(const string& dataset_name, int collapseType) {

    // Test merge function with items from dataset
    cout << endl << BOLDRED << "Test with dataset: " << dataset_name << " initial alpha = " << DEFAULT_ALPHA << " bin limit = " << DEFAULT_BIN_LIMIT << RESET << endl << endl;

    // init the sketch
    DDS_type *dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    DDS_type *dds2 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    // load dataset
    int n_element = getDatasetSize(dataset_name);

    // init array for all the elements
    double* dataset = nullptr;
    dataset = new (nothrow) double[n_element];
    if(!dataset){
        cout << "Not enough memory" << endl;
        DDS_Destroy(dds1);
        DDS_Destroy(dds2);
        return -2;
    }

    loadDataset(dataset_name, dataset, 0);

    cout << endl << BOLDRED << "Test with two sketches" << RESET << endl;

    // Insert items into two sketches
    insertFromDataset(dds1, 0, ceil(n_element/2), dataset, collapseType);
    insertFromDataset(dds2,(n_element/2), n_element, dataset, collapseType);
    merge(dds1,dds2,collapseType);

    testQuantile(dds1, dataset, n_element);
    deleteElements(dds1, dataset, n_element, collapseType);

    // reset dds1
    DDS_Destroy(dds1);
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    cout << BOLDRED << "Test with one sketch" << RESET << endl;

    // Test result with only one sketch
    insertFromDataset(dds1, 0, n_element, dataset, collapseType);
    testQuantile(dds1, dataset, n_element);
    deleteElements(dds1, dataset, n_element, collapseType);

    // deallocate the sketch data structure
    DDS_Destroy(dds1);
    DDS_Destroy(dds2);

    delete[] dataset, dataset = nullptr;

    return 0;
}

int testMergeFromTwoDataset(const string& dataset1_name, const string& dataset2_name, int collapseType) {

    // Test merge function with items from dataset
    cout << endl << BOLDRED << "Test with dataset: " << dataset1_name << " and " << dataset2_name << " initial alpha = " << DEFAULT_ALPHA << " bin limit = " << DEFAULT_BIN_LIMIT << RESET << endl << endl;

    // init the sketch
    DDS_type *dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    DDS_type *dds2 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    // load dataset
    int n_element1 = getDatasetSize(dataset1_name);
    int n_element2 = getDatasetSize(dataset2_name);

    // init array for all the elements
    double* dataset = nullptr;
    dataset = new (nothrow) double[n_element1+n_element2];
    if(!dataset){
        cout << "Not enough memory" << endl;
        DDS_Destroy(dds1);
        DDS_Destroy(dds2);
        return -2;
    }

    loadDataset(dataset1_name, dataset, 0);
    loadDataset(dataset2_name, dataset, n_element1 + 1);

    cout << endl << BOLDRED << "Test with two sketches" << RESET << endl;

    // Insert items into two sketches
    insertFromDataset(dds1, 0, n_element1, dataset, collapseType);
    insertFromDataset(dds2, n_element1, n_element1 + n_element2, dataset, collapseType);
    merge(dds1,dds2,collapseType);

    testQuantile(dds1, dataset, n_element1 + n_element2);
    deleteElements(dds1, dataset, n_element1 + n_element2, collapseType);

    // reset dds1
    DDS_Destroy(dds1);
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    cout << BOLDRED << "Test with one sketch" << RESET << endl;

    // Test result with only one sketch
    insertFromDataset(dds1, 0, n_element1 + n_element2, dataset, collapseType);
    testQuantile(dds1, dataset, n_element1 + n_element2);
    deleteElements(dds1, dataset, n_element1 + n_element2, collapseType);

    DDS_PrintCSV(dds1, "mergetwo.csv");

    // deallocate the sketch data structure
    DDS_Destroy(dds1);
    DDS_Destroy(dds2);

    delete[] dataset, dataset = nullptr;

    return 0;
}

int testQuantile(DDS_type *dds, double* stream, int n_element) {

    // Determine the quantile
    float q[] = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};
    int n_q = 11;

    cout << string(60, '-') << endl;
    cout << "|" << setw(10) << "quantile" << "|" << setw(15) << "estimate" << "|" << setw(15) << "real" << "|"  << setw(15)<< "error" << "|" << endl;
    cout << string(60, '-') << endl;

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

        cout << "|" << setw(10) << q[i]<< "|" << setw(15) << quantile << "|" << setw(15) << stream[idx-1] << "|"  << setw(15)<< error << "|" << endl;

    }

    cout << string(60, '-') << endl;

    return  0;
}

int deleteElements(DDS_type* dds, double* stream, int n_element, int collapseType) {

    // now check that delete works
    cout <<  endl << "Sketch size (number of bins) before delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl << endl;

    for (int i = 0; i < n_element; i++) {
        switch (collapseType) {
            case 1:
                DDS_DeleteCollapse(dds, stream[i]);
                break;
            case 2:
                DDS_DeleteCollapseLastBucket(dds, stream[i]);
                break;
            case 3:
                DDS_DeleteCollapseFirstBucket(dds, stream[i]);
        }
    }

    cout << "Sketch size (number of bins) after delete is equal to " << DDS_Size(dds) << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl << endl;

    return 0;
}

