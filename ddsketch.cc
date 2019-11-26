/********************************************************************
DDSketch

An algorithm for tracking quantiles in data streams

Charles Masson, Jee E. Rim, and Homin K. Lee. 2019. DDSketch: a fast and fully-mergeable quantile sketch with relative-error guarantees. Proc. VLDB Endow. 12, 12 (August 2019), 2195-2205. DOI: https://doi.org/10.14778/3352063.3352135

This implementation by
by Giuseppe Morleo
University of Salento, Italy

*********************************************************************/
/// \file

#include <fstream>
#include <iomanip>
#include "ddsketch.h"
#include "error.h"

high_resolution_clock::time_point t1, t2;

DDS_type *DDS_Init(int offset, int bin_limit, double alpha)
{

    // Initialize the sketch based on user-supplied parameters
    DDS_type *dds = nullptr;

    dds = new (nothrow) (DDS_type); // do not use the C malloc() function: it does not call C++ constructors
    if(!dds){
        fprintf(stdout,"Memory allocation of sketch data structure failed\n");
        return nullptr;
    }

    dds->offset = offset;
    dds->bin_limit = bin_limit;
    dds->alpha = alpha;
    dds->gamma = (1 + dds->alpha)/(1 - dds->alpha);
    dds->ln_gamma = log(dds->gamma);
    dds->bins = new map<int, int>();
    if(!dds->bins){
        printError(MEMORY_ERROR, __FUNCTION__);
        delete dds;
        return nullptr;
    }

    dds->n = 0;

    dds->max = numeric_limits<int>::min();
    dds->min = numeric_limits<int>::max();

    dds->min_value = pow(dds->gamma,pow(2,29));

    return dds;
}

void DDS_Destroy(DDS_type *dds)
{
    // get rid of a sketch and free up the space
    if (!dds) return;

    // deallocate the map
    if (dds->bins){
        delete dds->bins;
    }

    // now free the whole data structure
    delete dds;
}

int DDS_Size(DDS_type *dds)
{
    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // return the number of bins currently in the sketch
    return dds->bins->size();
}

int DDS_GetKey(DDS_type *dds, double item, int *error)
{

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // Given a value (item), returns the correspondning bucket index
    int key;

    if (item > 0) {
        key = int(ceil((log(item))/dds->ln_gamma)) + dds->offset;
    } else if (item < 0) {
        key = -int(ceil((log(-item))/dds->ln_gamma)) - dds->offset;
    } else if ( abs(item) < dds->min_value) {
        key = 0;
    }

    return key;

}

double DDS_GetRank(DDS_type *dds, int i, int *error)
{

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    double result;
    // Given a bucket index (i), this function returns the estimation of the rank x_q

    if (i > 0) {
        i -= dds->offset;
        result =  (2 * pow(dds->gamma, i))/(dds->gamma + 1);
    } else if (i < 0){
        i += dds->offset;
        result =  -(2 * pow(dds->gamma, -i))/(dds->gamma + 1);
    } else if ( i == 0) {
        result = dds->min_value;
    }

    return result;

}

double DDS_GetValue(DDS_type *dds, int i, int *error) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // if the bucket index is positive subtract the offset, otherwise add the offset
    // then, compute the bound by exponentiating gamma
    double bound;

    if ( i > 0) {
        i -= dds->offset;
        bound = pow(dds->gamma,i);
    } else {
        i += dds->offset;
        bound = -pow(dds->gamma,-i);
    }

    return bound;

}

int DDS_GetBounds(DDS_type *dds, int i, double* min, double* max) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // if the bucket index is positive subtract the offset, otherwise add the offset
    // then, compute the bound by exponentiating gamma

    if ( i > 0) {
        i -= dds->offset;
        (*max) = pow(dds->gamma,i);
        (*min) = pow(dds->gamma,i-1);
    } else {
        i += dds->offset;
        i = -i;
        *min = -pow(dds->gamma,i);
        *max = -pow(dds->gamma,i-1);
    }

    return 0;

}

int DDS_CollapseKey(DDS_type* dds, double i, int of, int *error){

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // if the bucket index is positive subtract the offset, otherwise add the offset
    // then, compute the new key

    if (i > 0) {
        i -= dds->offset;
        i = ceil((i+of)/2);
        i += dds->offset;
    } else if ( i< 0 ){
        i += dds->offset;
        i = -i;
        i = -ceil((i+of)/2);
        i -= dds->offset;
    } else if ( abs(i) < dds->min_value) {
        i = 0;
    }

    return int(i);
}


int DDS_AddCollapse(DDS_type *dds, double item)
{
    int *error = new int;
    *error = 0;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter

    int key = DDS_GetKey(dds, item, error);
    if ( *error < 0) {
        return -1;
    }

    (*(dds->bins))[key] += 1;
    dds->n += 1;

    int size = DDS_Size(dds);
    if ( size < 0 ) {
        return size;
    }

    while ( size > dds->bin_limit ) {
        // While the bin size is greater than bin_limit, we need to increase alpha and adapt all of the existing buckets to the new alpha value

        // collapse with gamma^2
        *error = DDS_Collapse(dds);
        if(*error < 0){
            return *error;
        }
    }

    return *error;

}

int DDS_AddCollapseLastBucket(DDS_type *dds, double item)
{
    int *error = new int;
    *error = 0;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter

    int key = DDS_GetKey(dds, item, error);
    if ( *error < 0) {
        return *error;
    }

    (*(dds->bins))[key] += 1;
    dds->n += 1;

    int size = DDS_Size(dds);
    if ( size < 0 ) {
        return size;
    }

    while ( size > dds->bin_limit ) {
        // If the bin size is greater than bin_limit, we need to increase alpha and adapt all of the existing buckets to the new alpha value

        // collapse the second last bucket into the last bucket
        *error = DDS_CollapseLastBucket(dds);
        if(*error < 0){
            return *error;
        }
    }

    return 0;

}

int DDS_AddCollapseFirstBucket(DDS_type *dds, double item)
{

    int *error = new int;
    *error = 0;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter

    int key = DDS_GetKey(dds, item, error);
    if ( *error < 0) {
        return *error;
    }

    (*(dds->bins))[key] += 1;
    dds->n += 1;

    int size = DDS_Size(dds);
    if ( size < 0 ) {
        return size;
    }

    while ( size > dds->bin_limit ) {
        // If the bin size is greater than bin_limit, we need to increase alpha and adapt all of the existing buckets to the new alpha value

        // collapse the second bucket into the first bucket
        *error = DDS_CollapseFirstBucket(dds);
        if(*error < 0){
            return *error;
        }
    }

    return *error;

}

int DDS_DeleteCollapse(DDS_type *dds, double item)
{
    int *error = new int;
    *error = 0;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // this function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1
    // otherwise it simply decrements by 1 the bucket's counter

    int key = DDS_GetKey(dds, item, error);
    if ( *error < 0) {
        return *error;
    }

    auto it = dds->bins->find(key);
    if (it != dds->bins->end()){

        // the bin associate to key actually exists
        // check its value: if it is 1 we erase the bin
        // otherwise we decrement by one the bin

        if(it->second == 1){
            dds->bins->erase(it);
            dds->n -= 1;
            //cout << "Deleted bin associated to item " << item << endl;
        }
        else{
            (*(dds->bins))[key] -= 1;
            dds->n -= 1;
            //cout << "Decremented bin associated to item " << item << endl;
        }


    }
    else{
        cout << item << ", " << key << " without offset " << DDS_RemoveOffset(dds, key, error) <<", \n";
        //dds->n -= 1;
        //cout << "There is no bin associated to item " << item << " with key " << key << endl;

    }

    return *error;

}

int DDS_DeleteCollapseLastBucket(DDS_type *dds, double item)
{
    int *error = new int;
    *error = 0;

    if(!dds){
        cout << "The sketch data structure is not valid" << endl;
        return -4;
    }

    // this function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1
    // otherwise it simply decrements by 1 the bucket's counter

    int key = DDS_GetKey(dds, item, error);
    if ( *error < 0) {
        return *error;
    }

    map<int, int>::iterator it;
    map<int, int>::reverse_iterator rit;

    if ( key >= dds->min && key <= dds->max ) {
        auto last = dds->bins->rbegin();
        key = last->first;
        it = dds->bins->find(key);
    } else {
        // the key is outside the [min,max] interval
        it = dds->bins->find(key);
    }

    if (it != dds->bins->end()){

        // the bin associate to key actually exists
        // check its value: if it is 1 we erase the bin
        // otherwise we decrement by one the bin

        if(it->second == 1){
            dds->bins->erase(it);
            dds->n -= 1;
            //cout << "Deleted bin associated to item " << item << endl;
        }
        else{
            (*(dds->bins))[key] -= 1;
            dds->n -= 1;
            //cout << "Decremented bin associated to item " << item << endl;
        }


    }
    else{
        //cout << item << ", " << DDS_RemoveOffset(dds, key) <<", \n";
        //dds->n -= 1;
        //cout << "There is no bin associated to item " << item << " with key " << key << endl;

    }

    return *error;

}

int DDS_DeleteCollapseFirstBucket(DDS_type *dds, double item)
{
    int *error = new int;
    *error = 0;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // this function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1
    // otherwise it simply decrements by 1 the bucket's counter

    int key = DDS_GetKey(dds, item, error);
    if ( *error < 0) {
        return *error;
    }

    map<int, int>::iterator it;

    if ( key >= dds->min && key <= dds->max ) {
        // the key within the [min,max] interval
        it = dds->bins->begin();
        key = it->first;
    } else {
        // the key is outside the [min,max] interval
        it = dds->bins->find(key);
    }

    if (it != dds->bins->end()){

        // the bin associate to key actually exists
        // check its value: if it is 1 we erase the bin
        // otherwise we decrement by one the bin

        if(it->second == 1){
            dds->bins->erase(it);
            dds->n -= 1;
            //cout << "Deleted bin associated to item " << item << endl;
        }
        else{
            (*(dds->bins))[key] -= 1;
            dds->n -= 1;
            //cout << "Decremented bin associated to item " << item << endl;
        }


    }
    else{
        //cout << item << ", " << DDS_RemoveOffset(dds, key) <<", \n";
        //dds->n -= 1;
        //cout << "There is no bin associated to item " << item << " with key " << key << endl;

    }

    return *error;

}

double DDS_GetQuantile(DDS_type *dds, float q, int *error)
{

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // If the value (q) is not in the [0,1] interval return NaN
    if (q < 0 || q > 1.01) {
        return numeric_limits<double>::quiet_NaN();
    }

    // We need to sum up the buckets until we find the bucket containing the q-quantile value x_q
    auto it = dds->bins->begin();
    int i = it->first;
    double count = it->second;
    double stop = q*double(dds->n - 1);

    while ( count <= stop) {

        ++it;
        i = it->first;
        count += it->second;

    }

    // Return the estimation x_q of bucket index i
    double rank = DDS_GetRank(dds, i, error);
    if ( *error < 0) {
        return *error;
    }

    return rank;

}

int DDS_MergeCollapse(DDS_type *dds1, DDS_type *dds2) {

    int *error = new int;
    *error = 0;

    if(!dds1 || !dds2){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    cout << "Size before merge sketch1 = " << DDS_Size(dds1) << " sketch2 = " << DDS_Size(dds2) << endl;

    startTheClock();

    // Check if the bins have the same alpha
    while (fabs(dds1->alpha - dds2->alpha) > 0.0000001){
        if (dds1->alpha < dds2->alpha) {
            *error = DDS_Collapse(dds1);
            if ( *error < 0 ) {
                return *error;
            }
        } else {
            *error = DDS_Collapse(dds2);
            if ( *error < 0 ) {
                return *error;
            }
        }
    }

    // Merge function merges the bins in dds1 with the bins of dds2
    // dds1 is the result of the merge operation
    for (auto received_bin : (*(dds2->bins))) {
        (*(dds1->bins))[received_bin.first] += received_bin.second;
        dds1->n += received_bin.second;
    }

    // Check if the new bin size is greater than bin limit
    int size = DDS_Size(dds1);
    if ( size < 0 ) {
        return size;
    }

    while ( size > dds1->bin_limit ) {
        // If the bin size is more then the bin limit, we need to collapse using gamma^2 instead of gamma
        *error = DDS_Collapse(dds1);
        if ( *error < 0 ) {
            return *error;
        }
    }

    double time = stopTheClock();

    cout << "Size after merge = " << DDS_Size(dds1) << " merge time = " << time << endl << endl;

    return *error;
}

int DDS_MergeCollapseLastBucket(DDS_type *dds1, DDS_type *dds2) {

    int *error = new int;
    *error = 0;

    if(!dds1 || !dds2){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    cout << "Size before merge sketch1 = " << DDS_Size(dds1) << " sketch2 = " << DDS_Size(dds2) << endl;

    startTheClock();

    // Check if the bins have the same alpha
    if (fabs(dds1->alpha - dds2->alpha) > 0.0000001){
        cout << "The two sketches cannot be merged, they have two different alphas" << endl;
        return MERGE_ERROR;
    }

    // Merge function merges the bins in dds1 with the bins of dds2
    // dds1 is the result of the merge operation
    for (auto received_bin : (*(dds2->bins))) {
        (*(dds1->bins))[received_bin.first] += received_bin.second;
        dds1->n += received_bin.second;
    }

    // Check if the new bin size is greater than bin limit
    int size = DDS_Size(dds1);
    if ( size < 0 ) {
        return size;
    }

    while (size > dds1->bin_limit){
        // If the bin size is more then the bin limit, we need to collapse the second last bucket in the last bucket
        *error = DDS_CollapseLastBucket(dds1);
        if ( *error < 0 ) {
            return *error;
        }
    }

    double time = stopTheClock();

    cout << "Size after merge = " << DDS_Size(dds1) << " merge time = " << time << endl << endl;

    return *error;
}

int DDS_MergeCollapseFirstBucket(DDS_type *dds1, DDS_type *dds2) {

    int *error = new int;
    *error = 0;

    if(!dds1 || !dds2){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    cout << "Size before merge sketch1 = " << DDS_Size(dds1) << " sketch2 = " << DDS_Size(dds2) << endl;

    startTheClock();

    // Check if the bins have the same alpha
    if (fabs(dds1->alpha - dds2->alpha) > 0.0000001){
        cout << "The two sketches cannot be merged, they have two different alphas" << endl;
        return MERGE_ERROR;
    }

    // Merge function merges the bins in dds1 with the bins of dds2
    // dds1 is the result of the merge operation
    for (auto received_bin : (*(dds2->bins))) {
        (*(dds1->bins))[received_bin.first] += received_bin.second;
        dds1->n += received_bin.second;
    }

    // Check if the new bin size is greater than bin limit
    int size = DDS_Size(dds1);
    if ( size < 0 ) {
        return size;
    }

    while (size > dds1->bin_limit){
        // If the bin size is more then the bin limit, we need to collapse the second bucket in the first bucket
        *error = DDS_CollapseFirstBucket(dds1);
        if ( *error < 0 ) {
            return *error;
        }
    }

    double time = stopTheClock();

    cout << "Size after merge = " << DDS_Size(dds1) << " merge time = " << time << endl << endl;

    return *error;
}

int DDS_CollapseLastBucket(DDS_type *dds) {

    int *error = new int;
    *error = 0;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // collapse the second bucket into the first bucket

    auto last = dds->bins->rbegin();
    auto second_last = --dds->bins->rbegin();

    if ( second_last->first < dds->min) {
        dds->min = (second_last->first);
    }
    if ( last->first > dds->max) {
        dds->max = last->first;
    }

    last->second += second_last->second;
    dds->bins->erase(second_last->first);

    return  0;
}

int DDS_CollapseFirstBucket(DDS_type *dds) {

    int *error = new int;
    *error = 0;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    // collapse the second last bucket into the last bucket

    auto first = dds->bins->begin();
    auto second = ++dds->bins->begin();

    if ( first->first < dds->min) {
        dds->min = (first->first);
    }
    if ( second->first > dds->max) {
        dds->max = (second->first);
    }

    first->second += second->second;
    dds->bins->erase(second->first);

    return  0;
}

int DDS_Collapse(DDS_type *dds) {

    int *error = new int;
    *error = 0;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    cout << "Size before collapse: " << DDS_Size(dds) << " sum: " << DDS_SumBins(dds) << " alpha: " << dds->alpha << " gamma: " << dds->gamma << endl;

    startTheClock();

    // determine new gamma and alpha values to be used
    dds->gamma = pow(dds->gamma, 2);
    dds->ln_gamma = log(dds->gamma);
    dds->alpha = (2 * dds->alpha) / (1 + pow(dds->alpha, 2));
    dds->min_value = pow(dds->gamma,pow(2,29));;

    // Create new bins map
    map<int, int> *new_bins = NULL;
    new_bins = new(nothrow) map<int, int>();
    if (!new_bins) {
        fprintf(stdout, "Memory allocation of a new sketch map failed\n");
        return -1;
    }

    // scan the buckets
    for (auto it = dds->bins->begin(); it != dds->bins->end(); ++it) {

        int key = it->first;
        // check if the bucket index is even
        if (key % 2 == 0) {
            int new_key = DDS_CollapseKey(dds, key, -1, error);
            if ( *error < 0 ) {
                return *error;
            }
            (*new_bins)[new_key] += it->second;
        } else {
            int new_key = DDS_CollapseKey(dds, key, +1, error);
            if ( *error < 0 ) {
                return *error;
            }
            (*new_bins)[new_key] += it->second;
        }

    }

    // Replace old bins map with new bins map
    dds->bins->swap(*new_bins);

    double time = stopTheClock();

    cout << "Size after collapse = " << DDS_Size(dds) << " sum: " << DDS_SumBins(dds) << " alpha: " << dds->alpha << " gamma: " << dds->gamma << " collapse time = " << time << endl << endl;

    new_bins->clear();
    delete new_bins;

    return 0;
}

int DDS_PrintCSV(DDS_type* dds, string name) {

    int *error = new int;
    *error = 0;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        *error = SKETCH_ERROR;
    }

    ofstream file;
    file.open(name);
    if (file.fail()) {
        printError(FILE_ERROR, __FUNCTION__);
        return FILE_ERROR;
    }

    file << fixed;
    file << setprecision(8);
    file << "key, count, max, min, length, \n";
    for (auto & bin : (*(dds->bins))) {

        auto* max = new double;
        auto* min = new double;
        DDS_GetBounds(dds, bin.first, min, max);

        file << DDS_RemoveOffset(dds, bin.first, error) <<", "<<bin.second<<", "<< *max<<", "<< *min <<", "<<max-min<<", \n";

    }

    file.close();

    return 0;
}

long DDS_SumBins(DDS_type *dds) {

    if(!dds){
        cout << "The sketch data structure is not valid" << endl;
        return -4;
    }

    long sum = 0;

    for (auto & bin : (*(dds->bins))) {
        sum += bin.second;
    }

    return sum;
}

int DDS_RemoveOffset(DDS_type* dds, int i, int *error) {

    if(!dds){
        cout << "The sketch data structure is not valid" << endl;
        return -4;
    }

    if (i > 0) {
        i -= dds->offset;
    } else {
        i += dds->offset;
    }

    return i;
}

int DDS_AddOffset(DDS_type* dds, int i){

    if(!dds){
        cout << "The sketch data structure is not valid" << endl;
        return -4;
    }

    if (i > 0) {
        i += dds->offset;
    } else {
        i -= dds->offset;
    }

    return i;
}

void startTheClock(){
    t1 = high_resolution_clock::now();
}

double stopTheClock() {
    t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    return time_span.count();
}
