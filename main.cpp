#include <iostream>
#include <math.h>
#include <limits>
#include <map>
#include <algorithm>
#include <random>

/*!
 * \mainpage Project $(project_name) Lorem ipsum dolor
 */
/// \file
using namespace std;

const int DEFAULT_OFFSET = 1073741824; //2^31/2
const int DEFAULT_BIN_LIMIT = 500;
const float DEFAULT_ALPHA = 0.01;

/*!
 * \brief               Given a value x, getKey returns the bucket index
 * @param x             The the input value
 * @param ln_gamma      The value of log(gamma)
 * @param offset        The offset to distinguish the values of negative x from positive x
 * @return              The index of the bucket containing the value x
 */
int getKey(double x, float ln_gamma, int offset);

/*!
 * \brief               Given a bucket's index, getRank returns the estimation of x_q
 * @param i             The index of the bucket
 * @param gamma         The basis of the logarithm related to alpha by the relation (1+alpha)/(1-alpha)
 * @param offset        The offset to distinguish the values of negative x from positive x
 * @return              The estimation of the x_q of that bucket
 */
double getRank(int i, float gamma, int offset);

/*!
 * \brief               The add function create a new bucket with index associated with the value x, or if that bucket already exists, it simply add 1 on the bucket's counter
 * @param x             The the input value to add to the sketch
 * @param bins          The hashmap containing the buckets
 * @param bin_limit     The maximum number of bucket
 * @param n             The global counter of all elements included in the sketch
 * @param alpha         The alpha-accuraxy level of q-quantile
 * @param gamma         The basis of the logarithm related to alpha by the relation (1+alpha)/(1-alpha)
 * @param ln_gamma      The value of log(gamma)
 * @param offset        The offset to distinguish the values of negative x from positive x
 */
void add(double x, map<int,int> &bins, int bin_limit, int &n, float &alpha, float &gamma, float &ln_gamma, int offset);

/*!
 * \brief               The expand function reduces the number of buckets, increasing the alpha, and then remapping each bucket according to the new range.
 * @param bins          The hashmap containing the buckets
 * @param alpha         The alpha-accuraxy level of q-quantile
 * @param gamma         The basis of the logarithm related to alpha by the relation (1+alpha)/(1-alpha)
 * @param ln_gamma      The value of log(gamma)
 * @param offset        The offset to distinguish the values of negative x from positive x
 */
void expand(map<int,int> &bins, float &alpha, float &gamma, float &ln_gamma, int offset);;

/*!
 * \brief               The getQuantile function returns the estimate of the desired q-quantile
 * @param q             The desired q-quantile
 * @param bins          The hashmap containing the buckets
 * @param n             The global counter of all elements included in the sketch
 * @param gamma         The basis of the logarithm related to alpha by the relation (1+alpha)/(1-alpha)
 * @param offset        The offset to distinguish the values of negative x from positive x
 * @return              The estimate of the desired q-quantile
 */
double getQuantile(float q, map<int,int> &bins, int n, float gamma, int offset);

/*!
 * \brief               The merge function merges the bins map with the received received_bins map
 * @param bins          The hashmap containing the buckets
 * @param received_bins The hashmap containing the sender's buckets
 * @param n             The global counter of all elements included in the sketch
 * @param bin_limit     The maximum number of bucket
 * @param alpha         The alpha-accuraxy level of q-quantile
 * @param gamma         The basis of the logarithm related to alpha by the relation (1+alpha)/(1-alpha)
 * @param ln_gamma      The value of log(gamma)
 * @param offset        The offset to distinguish the values of negative x from positive x
 */
void merge(map<int,int> &bins, const map<int,int> &received_bins, int &n, int bin_limit, float &alpha, float &gamma, float &ln_gamma, int offset);;

/// \brief  Main function
/// \param  argc An integer argument count of the command line arguments
/// \param  argv An argument vector of the command line arguments
/// \return an integer 0 upon exit success
int main() {

    // Parameters of the algorithm
    int offset = DEFAULT_OFFSET;
    int bin_limit = DEFAULT_BIN_LIMIT;
    float alpha = DEFAULT_ALPHA;

    float gamma = (1 + alpha)/(1-alpha);
    float ln_gamma = log(gamma);

    // Bins map
    map<int,int> bins;
    // Global counter of inserted elements
    int n = 0;

    // Test with normal distribution
    int n_element = 10000000;
    vector<double> stream;

    // generate pseudorandom number x according to the normal distribution
    // insert x into the stream vector and
    // add it to our DDSketch
    default_random_engine generator;
    normal_distribution<double> distribution(2,4);
    double x;
    for (int i=0; i<n_element;i++) {
        x = distribution(generator);
        stream.insert(stream.end(),x);
        add(x, bins,bin_limit, n, alpha, gamma, ln_gamma, offset);
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

    double quantile = getQuantile(q, bins, n, gamma, offset);
    double error = abs((quantile-stream[idx-1])/stream[idx-1]);

    cout << "Result for q=" << q << endl;
    cout << "Real: " << stream[idx-1] << " Estimation: " << quantile << " Error: " << error << endl;

    return 0;
}

int getKey(double x, float ln_gamma, int offset) {

    int key;

    if ( x > 0) {
        key = int(ceil((log(x))/ln_gamma)) + offset;
    } else if (x < 0) {
        key = -int(ceil((log(-x))/ln_gamma)) - offset;
    } else {
        key = 0;
    }

    return key;
}

double getRank(int i, float gamma, int offset) {

    if ( i > 0) {
        i = i - offset;
        return (2*pow(gamma,i))/(gamma+1);
    } else {
        i = i + offset;
        return -(2*pow(gamma,-i))/(gamma+1);
    }
}

void add(double x, map<int,int> &bins, int bin_limit, int &n, float &alpha, float &gamma, float &ln_gamma, int offset) {

    int key = getKey(x, ln_gamma, offset);
    bins[key] += 1;
    n += 1;

    if ( bins.size() > bin_limit ){
        // If the bin size is more than the bin limit, we need to increase alpha and adapt all the existing buckets with
        // the new alpha
        expand(bins, alpha, gamma, ln_gamma, offset);
    }
}

void expand(map<int,int> &bins, float &alpha, float &gamma, float &ln_gamma, int offset) {

    // We compute the new values of gamma and ln_gamma according the new alpha.
    alpha += 0.01;
    cout << "New alpha = " << alpha << endl;
    float new_gamma = (1 + alpha)/(1-alpha);
    float new_ln_gamma = log(new_gamma);

    double x;
    int key;

    // Create new bins map
    map<int,int> new_bins;

    for (auto & bin : bins) {
        x = getRank(bin.first, gamma, offset);
        key = getKey(x, new_ln_gamma, offset);
        new_bins[key] += bin.second;
    }

    // Replace old bins map with new bins map
    bins.swap(new_bins);
    new_bins.clear();

    // Replace gamma and ln_gamma with the new values according the new alpha
    gamma = new_gamma;
    ln_gamma = new_ln_gamma;
}

double getQuantile(float q, map<int,int> &bins, int n, float gamma, int offset) {

    // If q value is not in the [0,1] interval return NaN
    if (q < 0 || q > 1.01) {
        return numeric_limits<double>::quiet_NaN();
    }

    // We need to sum up the buckets until it finds the bucket containing the q-quantile value x_q
    auto it = bins.begin();
    int i = it->first;
    int count = it->second;
    while (count <= q*(n-1)) {
        ++it;
        i = it->first;
        count += it->second;
    }

    return getRank(i, gamma, offset);
}

void merge(map<int,int> &bins, const map<int,int> &received_bins, int &n, int bin_limit, float &alpha, float &gamma, float &ln_gamma, int offset) {

    for (auto received_bin : received_bins) {
        bins[received_bin.first] += received_bin.second;
        n += received_bin.second;
    }

    // After the merge operation we need to check if the new bin size is greater than bin limit
    if ( bins.size() > bin_limit ){
        // If the bin size is more then the bin limit, we need to increase alpha and adapt all the existing buckets with
        // the new alpha
        expand(bins, alpha, gamma, ln_gamma, offset);;
    }
}