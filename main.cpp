#include <iostream>
#include <math.h>
#include <limits>
#include <map>
#include <algorithm>
#include <random>

using namespace std;

const int DEFAULT_OFFSET = 1073741824; //2^31/2
const int DEFAULT_BIN_LIMIT = 500;
const float DEFAULT_ALPHA = 0.01;

int getKey(double x, float &ln_gamma, int &offset);
double getRank(int i, float &gamma, int &offset);
void add(double x, map<int,int> &bins, int &bin_limit, int &n, float &alpha, float &gamma, float &ln_gamma, int &offset);
void expand(map<int,int> &bins, float &alpha, float &gamma, float &ln_gamma, int &offset);
double getQuantile(float &q, map<int,int> &bins, int &n, float &gamma, int &offset);
void merge(map<int,int> &bins, const map<int,int> &received_bins, int &n, int &bin_limit, float &alpha, float &gamma, float &ln_gamma, int &offset);

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

    default_random_engine generator;
    normal_distribution<double> distribution(2,4);
    double x;
    for (int i=0; i<n_element;i++) {
        x = distribution(generator);
        stream.insert(stream.end(),x);
        add(x, bins,bin_limit, n, alpha, gamma, ln_gamma, offset);
    }

    sort(stream.begin(),stream.end());

    float q = 0.7;
    int idx = floor(1+q*(stream.size()-1));
    double quantile = getQuantile(q, bins, n, gamma, offset);
    double error = (quantile-stream[idx-1])/stream[idx-1];

    cout << "Result for q=" << q << endl;
    cout << "Real: " << stream[idx-1] << " Estimation: " << quantile << " Error: " << error << endl;

    return 0;
}

int getKey(double x, float &ln_gamma, int &offset) {
    // Given a value x, getKey returns the bucket index
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

double getRank(int i, float &gamma, int &offset) {
    // Given a bucket index, getRank returns the estimation of x_q
    if ( i > 0) {
        i = i - offset;
        return (2*pow(gamma,i))/(gamma+1);
    } else {
        i = i + offset;
        return -(2*pow(gamma,-i))/(gamma+1);
    }
}

void add(double x, map<int,int> &bins, int &bin_limit, int &n, float &alpha, float &gamma, float &ln_gamma, int &offset) {
    // The add function create a new bucket with index associated with the value x, or if that bucket already exists, it
    // simply add 1 on the bucket's counter
    int key = getKey(x, ln_gamma, offset);
    bins[key] += 1;
    n += 1;
    if ( bins.size() > bin_limit ){
        // If the bin size is more than the bin limit, we need to increase alpha and adapt all the existing buckets with
        // the new alpha
        expand(bins, alpha, gamma, ln_gamma, offset);
    }
}

void expand(map<int,int> &bins, float &alpha, float &gamma, float &ln_gamma, int &offset) {
    // In order to reduce the bucket's number, we need to increase the range on the bucket's index.
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

double getQuantile(float &q, map<int,int> &bins, int &n, float &gamma, int &offset) {
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
    // Return the estimation x_q of bucket index i
    return getRank(i, gamma, offset);
}

void merge(map<int,int> &bins, const map<int,int> &received_bins, int &n, int &bin_limit, float &alpha, float &gamma, float &ln_gamma, int &offset) {
    // Merge function merges our bins map with the received bins map
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