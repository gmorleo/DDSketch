/********************************************************************
 DDSketch
 
 An algorithm for tracking quantiles in data streams
 
 Charles Masson, Jee E. Rim, and Homin K. Lee. 2019. DDSketch: a fast and fully-mergeable quantile sketch with relative-error guarantees. Proc. VLDB Endow. 12, 12 (August 2019), 2195-2205. DOI: https://doi.org/10.14778/3352063.3352135

 This implementation by
 by Giuseppe Morleo
 University of Salento, Italy
 
 *********************************************************************/
/// \file

#include <iostream>
#include <math.h>
#include <limits>
#include <map>
#include <algorithm>
#include <random>

using namespace std;

typedef struct DDS_type{
    /// used to allow the skecth storing both positive, 0 and negative values
    int offset;
    /// maximum number of bins
    int bin_limit;
    /// this parameter defines alpha-accuracy of a q-quantile
    double alpha;
    /// this parameter is defined as (1 + alpha)/(1-alpha)
    double gamma;
    /// this is not a required parameter; it is defined as log(gamma)
    double ln_gamma;
    /// this map implements the bins of DDSketch
    map<int, int> *bins;
    /// this parameter keeps track of the number of items added to the sketch
    int n;

    /// collapse
    int min;
    int max;
} DDS_type;

/**
 * @brief               DDS costructor
 * @param offset        Used to allow the skecth storing both positive, 0 and negative values
 * @param bin_limit     The maximum number of bins
 * @param alpha         The alpha-accuraxy level of q-quantile
 * @return              an allocated DDSketch data structure
 */
extern DDS_type *DDS_Init(int offset, int bin_limit, double alpha);

/**
 * /brief               DDS destructor: deallocates a DDSketch data structure
 * @param dds           an allocated DDSketch data structure
 */
extern void DDS_Destroy(DDS_type *dds);
/**
 * \brief               Return the number of bins currently in the sketch
 * @param dds           The sketch
 * @return              The bins size, (number of bins)
 */
extern long DDS_Size(DDS_type *dds);

/**
 * \brief               Given a value (item), getKey returns the bucket index
 * @param dds           The sketch
 * @param item          The input value
 * @return              The index of the bucket containing the item
 */
extern int DDS_GetKey(DDS_type *dds, double item);

/**
 * \brief               Given a value (item), getKey returns the bucket index
 * @param dds           The sketch
 * @param item          The input value
 * @param ln_gamma      the log(gamma) to be used
 * @return              The index of the bucket containing the item
 */
extern int DDS_GetKey(DDS_type *dds, double item, float ln_gamma);

/**
 * \brief               Given a bucket index (i), this function returns the estimation of the rank x_q
 * @param dds           The sketch
 * @param i             The key of the bucket
 * @return              The estimate of the rank x_q
 */
extern double DDS_GetRank(DDS_type *dds, int i);

/**
 * \brief               This function creates a new bucket with index associated with the value (item), or if that bucket already exists, it simply add 1 to the bucket's counter
 * @param dds           The sketch
 * @param item          The the input value
 * @return              0 success, -1 error
 */
extern int DDS_AddCollapse(DDS_type *dds, double item);

extern int DDS_AddCollapseLastBucket(DDS_type *dds, double item);

/**
 * @brief               This function returns the bound associated with the key (i), (gamma^i)
 * @param dds           The sketch
 * @param i             The key of the bucket
 * @return
 */
double DDS_GetBound(DDS_type *dds, int i);

/**
 * @brief               This function returns the bound associated with the key (i), (gamma^i)
 * @param dds           The sketch
 * @param i             The key of the bucket
 * @param gamma         Gamma
 * @return
 */
double DDS_GetBound(DDS_type *dds, int i, float gamma);

/**
 * \brief               This function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1 otherwise it simply decrements by 1 the bucket's counter
 * @param dds           The sketch
 * @param item          The input value
 * @return              0 success
 */
extern int DDS_Delete(DDS_type *dds, double item);

extern int DDS_DeleteCollapseLastBucket(DDS_type *dds, double item);

/**
 * \brief               The function computes the estimate of the desired q-quantile (q)
 * @param dds           The sketch
 * @param q             The desired q-quantile
 * @return              The estimate of the desired q-quantile
 */
extern double DDS_GetQuantile(DDS_type *dds, float q);

/**
 * \brief               Merge function: merges the bins in dds1 with the bins of dds2; dds1 is the result of the merge operation
 * @param dds1          The sketch
 * @param dds2          The sketch
 */
extern void DDS_merge(DDS_type *dds1, DDS_type *dds2);

/**
 * @brief               This function computes the sum of the counters stored in the bins
 * @param dds           The sketch
 * @return              long
 */
extern long DDS_SumBins(DDS_type *dds);

/**
 * @brief               This function prints the bins map in a CSV file
 * @param dds           The sketch
 * @param name          File name
 * @return              0 success
 */
extern int DDS_PrintCSV(DDS_type* dds, string name);

/**
 * @brief               This function checks if a given item has a corresponding bucket
 * @param dds           The sketch
 * @param item          Input value
 * @return              0 success, -1 failure
 */
extern int DDS_CheckItem(DDS_type *dds, double item);

/**
 * @brief               This function collapses the last two buckets
 * @param dds           The sketch
 * @return
 */
extern int DDS_CollapseLastBucket(DDS_type *dds);

/**
 * @brief               The function collapses the old buckets in the new buckets based on the new range (range ^ 2)
 * @param dds           The sketch
 * @return              0 success, -1 failure
 */
extern int DDS_Collapse(DDS_type *dds);

/**
 * @brief               This function subtracts from the bucket index the offset used in the implementation to handle both positive and negative values
 * @param dds           The sketch
 * @param i             the bucket index
 * @return              the bucket index minus the offset
 */
extern int DDS_RemoveOffset(DDS_type* dds, int i);

/**
* @brief               This function adds to the bucket index the offset used in the implementation to handle both positive and negative values
* @param dds           The sketch
* @param i             the bucket index
* @return              the bucket index plus the offset
*/
extern int DDS_AddOffset(DDS_type* dds, int i);

/**
 *
 * @param dds           The sketch
 * @param i             the bucket index
 * @param of
 * @return              the bucket index
 */
extern int DDS_NewKey(DDS_type* dds,  double i, int of);