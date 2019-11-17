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



DDS_type *DDS_Init(int offset, int bin_limit, double alpha)
{

    // Initialize the sketch based on user-supplied parameters
    DDS_type *dds = NULL;

    dds = new (nothrow) (DDS_type); // do not use the C malloc() function: it does not call C++ constructors
    if(!dds){
        fprintf(stdout,"Memory allocation of sketch data structure failed\n");
        return NULL;
    }

    dds->offset = offset;
    dds->bin_limit = bin_limit;
    dds->alpha = alpha;
    dds->gamma = (1 + dds->alpha)/(1 - dds->alpha);
    dds->ln_gamma = log(dds->gamma);
    dds->bins = new map<int, int>();
    if(!dds->bins){
        fprintf(stdout,"Memory allocation of sketch map failed\n");
        delete dds;
        return NULL;
    }

    dds->n = 0;

    //remap
    dds->remap = new map<int, int>();
    if(!dds->remap){
        fprintf(stdout,"Memory allocation of sketch map failed\n");
        delete dds;
        return NULL;
    }

    dds->remapped = false;

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

    if (dds->remap){
        delete dds->remap;
    }

    // now free the whole data structure
    delete dds;
}

long DDS_Size(DDS_type *dds)
{
    // return the number of bins currently in the sketch
    return dds->bins->size();
}

int DDS_GetKey(DDS_type *dds, double item)
{
    // Given a value (item), returns the correspondning bucket index
    int key;

    if (item > 0) {
        key = int(ceil((log(item))/dds->ln_gamma)) + dds->offset;
    } else if (item < 0) {
        key = -int(ceil((log(-item))/dds->ln_gamma)) - dds->offset;
    } else {
        key = 0;
    }

    return key;

}

int DDS_GetKey(DDS_type *dds, double item, float ln_gamma)
{
    // Given a value (item), returns the correspondning bucket index
    int key;

    if (item > 0) {
        key = int(ceil((log(item))/ln_gamma)) + dds->offset;
    } else if (item < 0) {
        key = -int(ceil((log(-item))/ln_gamma)) - dds->offset;
    } else {
        key = 0;
    }

    return key;

}

double DDS_GetRank(DDS_type *dds, int i)
{

    // Given a bucket index (i), this function returns the estimation of the rank x_q

    if (i > 0) {
        i -= dds->offset;
        return (2 * pow(dds->gamma, i))/(dds->gamma + 1);
    } else {
        i += dds->offset;
        return -(2 * pow(dds->gamma, -i))/(dds->gamma + 1);
    }
}

double DDS_GetBound(DDS_type *dds, int i) {

    // This function returns the bound associated with the key (i), (gamma^i)

    if ( i > 0) {
        i -= dds->offset;
        return pow(dds->gamma,i);
    } else {
        i += dds->offset;
        return -pow(dds->gamma,-i);
    }

}

double DDS_GetBound(DDS_type *dds, int i, float gamma) {

    // This function returns the bound associated with the key (i), (gamma^i)

    if ( i > 0) {
        i -= dds->offset;
        return pow(gamma,i);
    } else {
        i += dds->offset;
        return -pow(gamma,-i);
    }

}

int DDS_Add(DDS_type *dds, double item)
{

    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter

    int key = DDS_GetKey(dds, item);
    (*(dds->bins))[key] += 1;
    dds->n += 1;

    if (DDS_Size(dds) > dds->bin_limit ){
        // If the bin size is greater than bin_limit, we need to increase alpha and adapt all of the existing buckets to the new alpha value

        // expand with alpha += 0.01
        //int return_value = DDS_expand(dds);

        // expand with alpha + = 0.01 with proportional redistribution
        //int return_value = DDS_expandProportional(dds);

        // collapse the last two bins
        //int return_value = DDS_Collapse(dds);

        // collapse with gamma^2
        int return_value = DDS_CollapsePlus(dds);
        if(return_value < 0){
            return -1;
        }
    }

    return 0;

}

int DDS_AddRemapped(DDS_type *dds, double item)
{

    // The Add function according the DDS_CollapseNeighbors

    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter

    int key = DDS_GetKey(dds, item);

    if ( dds->remapped && (dds->remap->find(key) != dds->remap->end())) {
        key = (*dds->remap)[key];
    }

    (*(dds->bins))[key] += 1;
    dds->n += 1;

    if (DDS_Size(dds) > dds->bin_limit ){
        // If the bin size is greater than bin_limit, we need to increase alpha and adapt all of the existing buckets to the new alpha value

        int return_value = DDS_CollapseNeighbors(dds);
        if(return_value < 0){
            return -1;
        }

    }

    return 0;

}

int DDS_Delete(DDS_type *dds, double item)
{

    // this function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1
    // otherwise it simply decrements by 1 the bucket's counter

    int key = DDS_GetKey(dds, item);

    map<int, int>::iterator it = dds->bins->find(key);
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
        return -1;
        //cout << item << ", " << DDS_RemoveOffset(dds, key) <<", \n";
        //dds->n -= 1;
        //cout << "There is no bin associated to item " << item << " with key " << key << endl;

    }

    return 0;

}

int DDS_DeleteCollapseNeighbors(DDS_type *dds, double item)
{

    // Delete function according the DDS_AddRemapped

    // this function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1
    // otherwise it simply decrements by 1 the bucket's counter

    int key = DDS_GetKey(dds, item);

    map<int, int>::iterator it;
    if ( dds->remap->find(key) != dds->remap->end()) {
        it = dds->bins->find(dds->remap->at(key));
    } else {
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
            it->second -= 1;
            //(*(dds->bins))[key] -= 1;
            dds->n -= 1;
            //cout << "Decremented bin associated to item " << item << endl;
        }


    }
    else{
        //dds->n -= 1;
        //cout << "There is no bin associated to item " << item << " with key " << key << endl;

    }

    return 0;

}

int DDS_expand(DDS_type *dds)
{
    // In order to reduce the bucket's number, we need to increase the range of the bucket's index.
    // We compute the new values of gamma and ln_gamma according the new alpha.
    dds->alpha += 0.01;
    cout << "New alpha = " << dds->alpha << endl;
    double new_gamma = (1 + dds->alpha)/(1 - dds->alpha);
    double new_ln_gamma = log(new_gamma);

    double item;
    int key;

    // Create new bins map
    map<int,int> *new_bins = NULL;
    new_bins = new (nothrow) map<int, int>();
    if(!new_bins){
        fprintf(stdout,"Memory allocation of a new sketch map failed\n");
        return -1;
    }

    for (auto & bin : (*(dds->bins))) {

        item = DDS_GetRank(dds, bin.first);
        key = DDS_GetKey(dds, item, new_ln_gamma);
        (*new_bins)[key] += bin.second;

    }

    // Replace old bins map with new bins map
    dds->bins->swap(*new_bins);
    new_bins->clear();
    delete new_bins;

    dds->gamma = new_gamma;
    dds->ln_gamma = new_ln_gamma;

    return 0;

}

int DDS_SumBins(DDS_type *dds) {

    // This function computes the sum of all counter of the bins

    int sum = 0;

    for (auto & bin : (*(dds->bins))) {
        sum += bin.second;
    }

    return sum;
}

int DDS_PrintCSV(DDS_type* dds, string name) {

    // This function prints the bins map in a CSV file

    ofstream file;
    file.open(name);
    if (file.fail()) {
        cout << "Error with file " << name << endl;
        return -1;
    }

    file << fixed;
    file << setprecision(8);
    file << "key, count, max, min, length, \n";

    for (auto & bin : (*(dds->bins))) {

        double max = DDS_GetBound(dds, bin.first);
        double min = DDS_GetBound(dds, (bin.first - 1));

        file << DDS_RemoveOffset(dds, bin.first) <<", "<<bin.second<<", "<<max<<", "<< min <<", "<<max-min<<", \n";

    }

    file.close();

    return 0;
}

int DDS_CheckAll(DDS_type *dds, double item) {

    // This function checks if all the elements in the stream have a corresponding bucket

    int key = DDS_GetKey(dds, item);

    auto it = dds->bins->find(key);
    if (it == dds->bins->end()){
        cout << "Not found key = " << key << endl;
        return -1;
    }

    return 0;
}

double DDS_GetQuantile(DDS_type *dds, float q)
{
    // If the value (q) is not in the [0,1] interval return NaN
    if (q < 0 || q > 1.01) {
        return numeric_limits<double>::quiet_NaN();
    }

    // We need to sum up the buckets until we find the bucket containing the q-quantile value x_q
    auto it = dds->bins->begin();
    int i = it->first;
    int count = it->second;

    while (count <= q*(dds->n - 1)) {

        ++it;
        i = it->first;
        count += it->second;

    }

    // Return the estimation x_q of bucket index i
    return DDS_GetRank(dds, i);

}

void DDS_merge(DDS_type *dds1, DDS_type *dds2)
{

    // Merge function merges the bins in dds1 with the bins of dds2
    // dds1 is the result of the merge operation
    for (auto received_bin : (*(dds2->bins))) {
        (*(dds1->bins))[received_bin.first] += received_bin.second;
        dds1->n += received_bin.second;
    }

    cout << "Size after merge = " << DDS_Size(dds1) << endl;

    // Check if the new bin size is greater than bin limit
    if (DDS_Size(dds1) > dds1->bin_limit){
        // If the bin size is more then the bin limit, we need to increase alpha and adapt all the existing buckets with
        // the new alpha
        DDS_expandProportional(dds1);
    }
}

int DDS_expandProportional(DDS_type *dds){

    // In order to reduce the bucket's number, we need to increase the range of the bucket's index.
    // We compute the new values of gamma and ln_gamma according the new alpha.
    dds->alpha += 0.01;
    double new_gamma = ((1 + dds->alpha)/(1 - dds->alpha));
    double new_ln_gamma = log(new_gamma);

    double item;
    int key;

    // Create new bins map
    map<int,int> *new_bins = NULL;
    new_bins = new (nothrow) map<int, int>();
    if(!new_bins){
        fprintf(stdout,"Memory allocation of a new sketch map failed\n");
        return -1;
    }

    for (auto & bin : (*dds->bins)) {

        item = DDS_GetRank(dds, bin.first);
        key = DDS_GetKey(dds, item, new_ln_gamma);

        double old_max = DDS_GetBound(dds, bin.first);
        double old_min = DDS_GetBound(dds, (bin.first - 1));

        double new_max = DDS_GetBound(dds, key, new_gamma);
        double new_min = DDS_GetBound(dds, (key - 1), new_gamma);

        double perc_next = 0;
        double perc_prec = 0;
        int current = 0;

        // If the old max is greater than the new max, a portion of elements must be distributed in the next bin
        if ( old_max > new_max) {

            perc_next = abs((old_max-new_max)/(old_max-old_min));

            // Saturation
            if ( perc_next > 1 ) {
                perc_next = 1;
            }

            (*new_bins)[key + 1] += floor(bin.second * perc_next);;

        }

        // If the old min is fewer than the new min, a portion of elements must be distributed in the precedent bin
        if ( old_min < new_min ) {

            perc_prec = abs((new_min-old_min)/(old_max-old_min));

            // Saturation
            if ( perc_prec > 1 ) {
                perc_prec = 1;
            }

            (*new_bins)[key - 1] += floor(bin.second * perc_prec);;

        }

        current = bin.second - floor(bin.second * perc_prec) - floor(bin.second * perc_next);

        if ( current != 0 ) {
            (*new_bins)[key] += current;
        }
    }

    // Replace old bins map with new bins map
    dds->bins->swap(*new_bins);
    new_bins->clear();
    delete new_bins;

    cout << "Size after expand = " << DDS_Size(dds) << endl;

    dds->gamma = new_gamma;
    dds->ln_gamma = new_ln_gamma;

    return 0;
}

int DDS_Collapse(DDS_type *dds) {

    // collapse the second last bucket into the last bucket

    auto last = dds->bins->rbegin();
    auto second_last = --dds->bins->rbegin();

    last->second += second_last->second;
    dds->bins->erase(second_last->first);

    return  0;
}

int DDS_CollapseNeighbors(DDS_type *dds){

    // It works for only one collapsing for now

    dds->remapped = true;

    // Create new bins map
    map<int,int> *new_bins = NULL;
    new_bins = new (nothrow) map<int, int>();
    if(!new_bins){
        fprintf(stdout,"Memory allocation of a new sketch map failed\n");
        return -1;
    }

    DDS_SumBins(dds);
    cout << "Size before expand = " << DDS_Size(dds) << endl;

    int prec_key = 0;
    for (auto & bin: (*dds->bins)) {

        int key = bin.first;
        if ( prec_key == key-1 ) {

            (*dds->remap)[key-1] = key-1;
            (*dds->remap)[key] = key-1;
            (*new_bins)[key-1] += bin.second;
            prec_key = 0;

        } else {

            (*new_bins)[key] += bin.second;
            prec_key = key;

        }
    }

    // Replace old bins map with new bins map
    dds->bins->swap(*new_bins);
    new_bins->clear();
    delete new_bins;

    cout << "Size after expand = " << DDS_Size(dds) << endl;

    return 0;

}

int DDS_RemoveOffset(DDS_type* dds, int i) {

    // This function removes the offset to the key (i)

    if (i > 0) {
        i -= dds->offset;
    } else {
        i += dds->offset;
    }
    return i;
}

int DDS_AddOffset(DDS_type* dds, int i){

    // This function adds the offset to the key (i)

    if (i > 0) {
        i += dds->offset;
    } else {
        i -= dds->offset;
    }
    return i;
}

int DDS_NewKey(DDS_type* dds,  double i, int of){


    if (i > 0) {
        i -= dds->offset;
        i = ceil((i+of)/2);
        i += dds->offset;
    } else {
        i += dds->offset;
        i = floor((i+of)/2);
        i -= dds->offset;
    }
    return int(i);
}

int DDS_CollapsePlus(DDS_type *dds) {

    cout << "Size before collapse: " << DDS_Size(dds) << " sum: " << DDS_SumBins(dds) << endl;

    dds->gamma = pow(dds->gamma, 2);
    dds->ln_gamma = log(dds->gamma);
    dds->alpha = (2 * dds->alpha) / (1 + pow(dds->alpha, 2));

    // Create new bins map
    map<int, int> *new_bins = NULL;
    new_bins = new(nothrow) map<int, int>();
    if (!new_bins) {
        fprintf(stdout, "Memory allocation of a new sketch map failed\n");
        return -1;
    }

    for (auto it = dds->bins->begin(); it != dds->bins->end(); ++it) {

        int key = it->first;
        if ( key > 0 ) {
            if ( key%2 == 0 ) {
                int new_key = DDS_NewKey(dds, key, -1);
                (*new_bins)[new_key] += it->second;
            } else {
                // check if the next bucket exist
                int new_key = DDS_NewKey(dds, key, +1);
                auto next_bin = next(it, 1);

                if ( next_bin->first == key+1) {
                    (*new_bins)[new_key] += it->second + next_bin->second;
                    // we have to skip the next bucket because we have already considered it in this interaction
                    ++it;
                    if ( it == dds->bins->end() ) {
                        break;
                    }
                } else {
                    (*new_bins)[new_key] += it->second;
                }
            }
        } else {
            if ( key%2 == 0 ) {
                // check if the next bucket exist
                int new_key = DDS_NewKey(dds, key, +1);
                auto next_bin = next(it, 1);

                if ( next_bin->first == key+1) {
                    (*new_bins)[new_key] += it->second + next_bin->second;
                    // we have to skip the next bucket because we have already considered it in this interaction
                    ++it;
                    if ( it == dds->bins->end() ) {
                        break;
                    }
                } else {
                    (*new_bins)[new_key] += it->second;
                }
            } else {
                int new_key = DDS_NewKey(dds, key, -1);
                (*new_bins)[new_key] += it->second;
            }
        }
    }

    // Replace old bins map with new bins map
    dds->bins->swap(*new_bins);
    new_bins->clear();
    delete new_bins;

    cout << "Size after collapse = " << DDS_Size(dds) << " sum: " << DDS_SumBins(dds) << " alpha: " << dds->alpha << " gamma: " << dds->gamma << endl;

    return 0;
}