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
#include "error.h"

#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define RESET   "\033[0m"

using namespace std;

const int DEFAULT_OFFSET = 1073741824; //2^31/2
const int DEFAULT_BIN_LIMIT = 500;
const double DEFAULT_ALPHA = 0.008;

int getErrorBound(DDS_type *dds, int collapseType);

/**
 * \brief                   This function computes the dimension of the dataset
 * @param name_file         Name of dataset
 * @return                  Return the number of element in the dataset(row)
 */
int getDatasetSize(const string &name_file, int &rows);

/**
 * \brief                   This function loads the dataset into an array
 * @param name_file         Name of dataset
 * @param dataset           Array
 * @return                  An array containing the whole dataset
 */
int loadDataset(const string &name_file, double *dataset);

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
template <class Distribution>
int printDataset(const string& name, int n_element, Distribution distribution);

/**
 * @brief                   This function inserts (n_element) to the sketch according to the selected distribution and collapse type
 * @param dds               The sketch
 * @param stream            Vector that contains all the real values inserted
 * @param stream_start      Initial index from which insert the elements in the vector
 * @param n_element         Number of element
 * @param distribution      Distribution
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
template <class Distribution>
int insertRandom(DDS_type *dds, double* stream, int n_element, Distribution distribution, int collapseType);

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
 * @param distribution      Distribution
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
template <class Distribution>
int testWithRandomValue(int n_element, Distribution distribution, int collapseType);

/**
 * @brief                   This function tests the merge of two sketches, entering (n_element1) in the first sketch and
 *                          (n_element2) in the second sketch based on the selected distributions and using the selected type of collapsing
 * @param n_element1        Number of element first distribution
 * @param n_element2        Number of element second distribution
 * @param distribution1     First distribution
 * @param distribution2     Second distribution
 * @param collapseType      Collapse type: 1 gamma^2, 2 last bucket, 3 first bucket
 * @return                  0 success -1 error
 */
template <class Distribution>
int testMergeWithRandomValue(int n_element1, int n_element2, Distribution distribution1, Distribution distribution2, int collapseType);

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

    /// Init the distribution
    default_random_engine generator;
    normal_distribution<double> normal(2,3);
    normal_distribution<double> normal2(10,3);
    exponential_distribution<double> exponential(17);
    uniform_real_distribution<double> uniform_real(-5,5);
    uniform_real_distribution<double> uniform_real2(0,10);
    gamma_distribution<double> gamma(2,2);

    /// number of element
    int n_element = pow(10,7);

    /// Test with random value
    /// Collapse type: 1 Collapse with gamma^2, 2 Collapse with last buckets, 3 Collapse with first buckets

    //testWithRandomValue(n_element, uniform_real, 3);
    //testWithRandomValue(n_element, uniform_real, 3);
    //testWithRandomValue(n_element, normal, 1);
    //testWithRandomValue(n_element, normal, 3);
    //testMergeWithRandomValue(n_element,n_element,uniform_real,uniform_real2,1);
    //testMergeWithRandomValue(n_element,n_element,uniform_real,uniform_real2,3);

    /// Test merge function
    /// Collapse type: 1 Collapse with gamma^2, 2 Collapse with last buckets, 3 Collapse with first bucket

    //testMergeFromDataset("normal.csv", 1);
    testMergeFromTwoDataset("normal.csv", "uniform_real.csv", 1);

    /// Print dataset on a file
    //printDataset("uniform_real.csv", 10000000, uniform_real);

    return 0;
}

template <class Distribution>
int printDataset(const string& name, int n_element, Distribution distribution) {

    // open file for output
    ofstream file;
    file.open(name);
    if (file.fail()) {
        printError(FILE_ERROR, __FUNCTION__);
        return FILE_ERROR;
    }

    // Init the distribution
    default_random_engine generator;

    double item;

    for (int i = 0; i < n_element; i++) {

        // generate pseudorandom number (item) according to the select distribution
        item = distribution(generator);

        // insert number (item) to the file
        file << item << ", \n";

    }

    file.close();

    return SUCCESS;
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

template <class Distribution>
int insertRandom(DDS_type *dds, double* stream, int n_element, Distribution distribution, int collapseType) {

    int returnValue = -1;

    // Init the distribution
    default_random_engine generator;

    double item;

    // generate pseudorandom number (item) according to the selected distribution
    // insert item into the stream vector and
    // add it to our DDSketch
    for (int i = 0; i < n_element; i++) {

        item = distribution(generator);

        stream[i] = item;

        switch (collapseType) {
            case 1:
                returnValue = DDS_AddCollapse(dds, item);
                if (returnValue) {
                    return  returnValue;
                }
                break;
            case 2:
                returnValue = DDS_AddCollapseLastBucket(dds, item);
                if (returnValue) {
                    return  returnValue;
                }
                break;
            case 3:
                returnValue = DDS_AddCollapseFirstBucket(dds, item);
                if (returnValue) {
                    return  returnValue;
                }
                break;
            default:
                printError(UNKNOWN_COLLAPSE_TYPE, __FUNCTION__);
                return UNKNOWN_COLLAPSE_TYPE;
        }

    }

    return returnValue;
}

int insertFromDataset(DDS_type* dds, int start, int stop, double* dataset, int collapseType) {

    int returnValue = -1;

    // insert item into  our DDSketch
    for (int i = start; i < stop; i++) {
        switch (collapseType) {
            case 1:
                returnValue = DDS_AddCollapse(dds, dataset[i]);;
                if (returnValue) {
                    return returnValue;
                }
                break;
            case 2:
                returnValue = DDS_AddCollapseLastBucket(dds, dataset[i]);
                if (returnValue) {
                    return returnValue;
                }
                break;
            case 3:
                returnValue = DDS_AddCollapseFirstBucket(dds, dataset[i]);
                if (returnValue) {
                    return returnValue;
                }
                break;
            default:
                printError(UNKNOWN_COLLAPSE_TYPE, __FUNCTION__);
                return UNKNOWN_COLLAPSE_TYPE;
        }
    }


    return returnValue;
}

int merge(DDS_type* dds1, DDS_type* dds2, int collapseType) {

    int returnValue = -1;

    // selects the merge function based on the selected collapse method
    switch (collapseType) {
        case 1:
            returnValue = DDS_MergeCollapse(dds1, dds2);
            if (returnValue) {
                return returnValue;
            }
            break;
        case 2:
            returnValue = DDS_MergeCollapseLastBucket(dds1, dds2);
            if (returnValue) {
                return returnValue;
            }
            break;
        case 3:
            returnValue = DDS_MergeCollapseFirstBucket(dds1, dds2);
            if (returnValue) {
                return returnValue;
            }
            break;
        default:
            printError(UNKNOWN_COLLAPSE_TYPE, __FUNCTION__);
            return UNKNOWN_COLLAPSE_TYPE;
    }

    return returnValue;
}

template <class Distribution>
int testWithRandomValue(int n_element, Distribution distribution, int collapseType) {

    /*** Test sketch with random values ***/

    cout << endl << BOLDRED << "Test with distribution: " << " initial alpha = " << DEFAULT_ALPHA << " bin limit = " << DEFAULT_BIN_LIMIT << RESET << endl;

    int returnValue = -1;
    DDS_type* dds1 = nullptr;
    double* stream = nullptr;

    /*** Init sketch ***/
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    if (!dds1) {
        goto ON_EXIT;
    }

    /*** Init array ***/
    stream = new (nothrow) double[n_element];
    if(!stream){
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    /*** Test with one sketch ***/

    // insert items into the sketch
    returnValue = insertRandom(dds1, stream, n_element, distribution, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // compute the range of wrong quantiles
    returnValue = getErrorBound(dds1, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test quantile
    returnValue = testQuantile(dds1, stream, n_element);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test delete
    returnValue = deleteElements(dds1, stream, n_element, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    ON_EXIT:

    if (dds1 != nullptr) {
        DDS_Destroy(dds1);
    }
    if (stream != nullptr) {
        delete[] stream, stream = nullptr;
    }

    return returnValue;
}

template <class Distribution>
int testMergeWithRandomValue(int n_element1, int n_element2, Distribution distribution1, Distribution distribution2, int collapseType) {

    /*** Test merge function with random values ***/

    cout << endl << BOLDRED << "Test with distribution: " << " initial alpha = " << DEFAULT_ALPHA << " bin limit = " << DEFAULT_BIN_LIMIT << RESET << endl;

    int returnValue = -1;
    DDS_type *dds1 = nullptr;
    DDS_type *dds2 = nullptr;
    double *stream1 = nullptr;
    double *stream2 = nullptr;
    double *stream = nullptr;

    /*** Init sketch ***/
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    if (!dds1) {
        goto ON_EXIT;
    }

    dds2 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    if (!dds2) {
        goto ON_EXIT;
    }

    /*** Init array ***/
    stream1 = new (nothrow) double[n_element1];
    if(!stream1){
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    stream2 = new (nothrow) double[n_element2];
    if(!stream2){
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    /*** Test with 2 sketches ***/

    cout << endl << BOLDRED << "Test with two sketches" << RESET << endl;

    // insert items into the first sketch
    returnValue = insertRandom(dds1, stream1, n_element1, distribution1, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // insert items into the second sketch
    returnValue = insertRandom(dds2, stream2, n_element2, distribution2, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Merge ***/

    returnValue = merge(dds1, dds2, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    stream = new (nothrow) double[n_element1+n_element2];
    if(!stream){
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    // Merge stream

    try {
        move( stream1, stream1+n_element1, stream );
        move( stream2, stream2+n_element2, stream + n_element1 );
    } catch (int e) {
        printError(COPY_ERROR, __FUNCTION__);
        returnValue = e;
        goto ON_EXIT;
    }

    // compute the range of wrong quantiles
    returnValue = getErrorBound(dds1, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test quantile
    returnValue = testQuantile(dds1, stream, n_element1 + n_element2);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test delete
    returnValue = deleteElements(dds1, stream, n_element1 + n_element2, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Test with one sketch ***/

    cout << BOLDRED << "Test with one sketch" << RESET << endl;

    DDS_Destroy(dds1);
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    if (!dds1) {
        goto ON_EXIT;
    }

    // insert items into the sketch
    returnValue = insertFromDataset(dds1, 0, n_element1 + n_element2, stream, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // compute the range of wrong quantiles
    returnValue = getErrorBound(dds1, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test quantile
    returnValue = testQuantile(dds1, stream, n_element1 + n_element2);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test delete
    returnValue = deleteElements(dds1, stream, n_element1 + n_element2, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Deallocate the sketch data structure ***/

    ON_EXIT:

    if (dds1 != nullptr) {
        DDS_Destroy(dds1);
    }

    if (dds2 != nullptr) {
        DDS_Destroy(dds2);
    }

    if (stream != nullptr) {
        delete[] stream, stream = nullptr;
    }

    if (stream1 != nullptr) {
        delete[] stream1, stream1 = nullptr;
    }

    if (stream2 != nullptr) {
        delete[] stream2, stream2 = nullptr;
    }

    return returnValue;
}

int testMergeFromDataset(const string& dataset_name, int collapseType) {

    /*** Test merge function with items from dataset ***/

    cout << endl << BOLDRED << "Test with dataset: " << dataset_name << " initial alpha = " << DEFAULT_ALPHA << " bin limit = " << DEFAULT_BIN_LIMIT << RESET << endl << endl;

    int returnValue = -1;
    DDS_type *dds1 = nullptr;
    DDS_type *dds2 = nullptr;
    double *dataset = nullptr;

    int n_element;

    /*** Init sketch ***/
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    if (!dds1) {
        goto ON_EXIT;
    }

    dds2 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    if (!dds2) {
        goto ON_EXIT;
    }

    /*** Init array ***/
    returnValue = getDatasetSize(dataset_name, n_element);
    if (returnValue) {
        goto ON_EXIT;
    }

    dataset = new (nothrow) double[n_element];
    if(!dataset){
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    /*** Load dataset ***/
    returnValue = loadDataset(dataset_name, dataset);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Test with 2 sketches ***/

    cout << endl << BOLDRED << "Test with two sketches" << RESET << endl;

    // insert items into the first sketch
    returnValue = insertFromDataset(dds1, 0, ceil(n_element/2), dataset, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // insert items into the second sketch
    returnValue = insertFromDataset(dds2,ceil(n_element/2), n_element, dataset, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Merge ***/

    returnValue = merge(dds1,dds2,collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // compute the range of wrong quantiles
    returnValue = getErrorBound(dds1, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test quantile
    returnValue = testQuantile(dds1, dataset, n_element);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test delete
    returnValue = deleteElements(dds1, dataset, n_element, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Test with one sketch ***/
    cout << BOLDRED << "Test with one sketch" << RESET << endl;

    DDS_Destroy(dds1);
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    if (!dds1) {
        goto ON_EXIT;
    }

    // insert items into the sketch
    returnValue = insertFromDataset(dds1, 0, n_element, dataset, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // compute the range of wrong quantiles
    returnValue = getErrorBound(dds1, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test quantile
    returnValue = testQuantile(dds1, dataset, n_element);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test delete
    returnValue = deleteElements(dds1, dataset, n_element, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Deallocate the sketch data structure ***/

    ON_EXIT:

    if (dds1 != nullptr) {
        DDS_Destroy(dds1);
    }

    if (dds2 != nullptr) {
        DDS_Destroy(dds2);
    }

    if (dataset != nullptr) {
        delete[] dataset, dataset = nullptr;
    }

    return returnValue;
}

int testMergeFromTwoDataset(const string& dataset1_name, const string& dataset2_name, int collapseType) {

    /*** Test merge function with items from dataset ***/

    cout << endl << BOLDRED << "Test with dataset: " << dataset1_name << " and " << dataset2_name << " initial alpha = " << DEFAULT_ALPHA << " bin limit = " << DEFAULT_BIN_LIMIT << RESET << endl << endl;

    int returnValue = -1;
    DDS_type *dds1 = nullptr;
    DDS_type *dds2 = nullptr;
    double *dataset = nullptr;
    double *dataset1 = nullptr;
    double *dataset2 = nullptr;

    int n_element1;
    int n_element2;

    /*** Init sketch ***/
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    if (!dds1) {
        goto ON_EXIT;
    }

    dds2 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);
    if (!dds2) {
        goto ON_EXIT;
    }

    /*** Init array ***/
    returnValue = getDatasetSize(dataset1_name, n_element1);
    if (returnValue) {
        goto ON_EXIT;
    }

    returnValue = getDatasetSize(dataset2_name, n_element2);
    if (returnValue) {
        goto ON_EXIT;
    }

    dataset1 = new (nothrow) double[n_element1];
    if(!dataset1){
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    dataset2 = new (nothrow) double[n_element2];
    if(!dataset2){
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    /*** Load dataset ***/
    returnValue = loadDataset(dataset1_name, dataset1);
    if (returnValue) {
        goto ON_EXIT;
    }

    returnValue = loadDataset(dataset2_name, dataset2);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Test with two sketches ***/

    cout << endl << BOLDRED << "Test with two sketches" << RESET << endl;

    // insert items into the first sketch
    returnValue = insertFromDataset(dds1, 0, n_element1, dataset1, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // insert items into the second sketch
    returnValue = insertFromDataset(dds2, 0, n_element2, dataset2, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Merge ***/

    returnValue = merge(dds1,dds2,collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    dataset = new (nothrow) double[n_element1+n_element2];
    if(!dataset){
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    // Merge dataset

    try {
        move( dataset1, dataset1+n_element1, dataset );
        move( dataset2, dataset2+n_element2, dataset + n_element1 );
    } catch (int e) {
        printError(COPY_ERROR, __FUNCTION__);
        returnValue = e;
        goto ON_EXIT;
    }

    // compute the range of wrong quantiles
    returnValue = getErrorBound(dds1, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test quantile
    returnValue = testQuantile(dds1, dataset, n_element1 + n_element2);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test delete
    returnValue = deleteElements(dds1, dataset, n_element1 + n_element2, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Test with one sketch ***/

    cout << BOLDRED << "Test with one sketch" << RESET << endl;

    DDS_Destroy(dds1);
    dds1 = DDS_Init(DEFAULT_OFFSET, DEFAULT_BIN_LIMIT, DEFAULT_ALPHA);

    // insert items into the sketch
    returnValue = insertFromDataset(dds1, 0, n_element1 + n_element2, dataset, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // compute the range of wrong quantiles
    returnValue = getErrorBound(dds1, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test quantile
    returnValue = testQuantile(dds1, dataset, n_element1 + n_element2);
    if (returnValue) {
        goto ON_EXIT;
    }

    // test delete
    returnValue = deleteElements(dds1, dataset, n_element1 + n_element2, collapseType);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Deallocate the sketch data structure ***/

    ON_EXIT:

    if (dds1 != nullptr) {
        DDS_Destroy(dds1);
    }

    if (dds2 != nullptr) {
        DDS_Destroy(dds2);
    }

    if (dataset != nullptr) {
        delete[] dataset, dataset = nullptr;
    }

    if (dataset1 != nullptr) {
        delete[] dataset1, dataset1 = nullptr;
    }

    if (dataset2 != nullptr) {
        delete[] dataset2, dataset2 = nullptr;
    }

    return returnValue;
}

int testQuantile(DDS_type *dds, double* stream, int n_element) {

    int returnValue = -1;

    // Determine the quantile
    float q[] = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};
    int n_q = 11;

    cout << string(60, '-') << endl;
    cout << "|" << setw(10) << "quantile" << "|" << setw(15) << "estimate" << "|" << setw(15) << "real" << "|"  << setw(15)<< "error" << "|" << endl;
    cout << string(60, '-') << endl;

    for ( int i = 0; i < n_q; i++ ) {

        int idx = floor(1+q[i]*double(n_element-1));
        // determine the correct answer
        // i.e., the number stored at the index (idx-1) in the sorted permutation of the vector
        // note that we are not sorting the vector, we are using the quickselect() algorithm
        // which in C++ is available as std::nth_element
        nth_element(stream, stream + (idx-1), stream +  n_element);
        double quantile;
        returnValue = DDS_GetQuantile(dds, q[i], quantile);
        if (returnValue < 0 ) {
            return returnValue;
        }

        double error = abs((quantile-stream[idx-1])/stream[idx-1]);

        cout << "|" << setw(10) << q[i]<< "|" << setw(15) << quantile << "|" << setw(15) << stream[idx-1] << "|"  << setw(15)<< error << "|" << endl;

    }

    cout << string(60, '-') << endl;

    return  returnValue;
}

int deleteElements(DDS_type* dds, double* stream, int n_element, int collapseType) {

    int returnValue = -1;

    int size;
    returnValue = DDS_Size(dds, size);
    if (returnValue) {
        return returnValue;
    }

    // now check that delete works
    cout <<  endl << "Sketch size (number of bins) before delete is equal to " << size << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl << endl;

    for (int i = 0; i < n_element; i++) {
        switch (collapseType) {
            case 1:
                returnValue = DDS_DeleteCollapse(dds, stream[i]);
                if (returnValue) {
                    return returnValue;
                }
                break;
            case 2:
                returnValue = DDS_DeleteCollapseLastBucket(dds, stream[i]);
                if (returnValue) {
                    return returnValue;
                }
                break;
            case 3:
                returnValue = DDS_DeleteCollapseFirstBucket(dds, stream[i]);
                if (returnValue) {
                    return returnValue;
                }
                break;
            default:
                printError(UNKNOWN_COLLAPSE_TYPE, __FUNCTION__);
                return UNKNOWN_COLLAPSE_TYPE;
        }
    }

    returnValue = DDS_Size(dds, size);
    if (returnValue < 0 ) {
        return returnValue;
    }

    cout << "Sketch size (number of bins) after delete is equal to " << size << endl;
    cout << "Number of items in the sketch is equal to " << dds->n << endl << endl;

    return returnValue;
}


int getDatasetSize(const string &name_file, int &rows) {

    rows = 0;

    ifstream inputFile(name_file);
    string line;
    if(inputFile.is_open()){
        while (getline(inputFile, line))
            rows++;
    } else {
        printError(FILE_ERROR, __FUNCTION__);
        return FILE_ERROR;
    }

    return SUCCESS;
}

int loadDataset(const string &name_file, double *dataset) {

    if (!dataset) {
        printError(NULL_POINTER_ERROR, __FUNCTION__);
        return NULL_POINTER_ERROR;
    }

    ifstream inputFile(name_file);
    string line;
    int row = 0;

    if(inputFile.is_open()){
        while (getline(inputFile, line)) {

            try {
                dataset[row] = stod(line);
            }
            catch (int e) {
                printError(FILE_ERROR, __FUNCTION__);
                inputFile.close();
                return FILE_ERROR;
            }

            row++;
        }
    } else {
        printError(FILE_ERROR, __FUNCTION__);
        return FILE_ERROR;
    }

    inputFile.close();

    return SUCCESS;
}

int getErrorBound(DDS_type *dds, int collapseType) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    int returnValue = -1;

    double count;
    double result;
    double min, max;

    switch (collapseType) {
        case 1:
            cout << endl << "The estimates are always correct" << endl;
            returnValue = 0;
            break;
        case 2:
            count = dds->bins->rbegin()->second;
            result = count / dds->n;
            cout << endl << "The estimates are wrong with quantiles between " <<  (1-result) << " and " << 1 << endl;
            returnValue = DDS_GetBounds(dds,dds->min, dds->max, min, max);
            if (returnValue) {
                break;
            }
            cout << "The last bucket contains value between " << min << " and " << max << ", last bucket counter = " << count << endl;
            break;
        case 3:
            count = dds->bins->begin()->second;
            result = count / dds->n;
            cout << endl << "The estimates are wrong with quantiles between " <<  0 << " and " << result << endl;
            returnValue = DDS_GetBounds(dds,dds->min, dds->max, min, max);
            if (returnValue) {
                break;
            }
            cout << "The first bucket contains value between " << min << " and " << max << ", first bucket counter = " << count << endl;
            break;
        default:
            printError(UNKNOWN_COLLAPSE_TYPE, __FUNCTION__);
            returnValue = UNKNOWN_COLLAPSE_TYPE;
    }

    return returnValue;
}

