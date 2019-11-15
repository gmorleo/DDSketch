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
    float alpha;
    /// this parameter is defined as (1 + alpha)/(1-alpha)
    float gamma;
    /// this is not a required parameter; it is defined as log(gamma)
    float ln_gamma;
    /// this map implements the bins of DDSketch
    map<int, int> *bins;
    /// this parameter keeps track of the number of items added to the sketch
    int n;

    /// remap
    map<int, int> *remap;
    bool remapped;

} DDS_type;

/**
 * @brief               DDS costructor
 * @param offset        Used to allow the skecth storing both positive, 0 and negative values
 * @param bin_limit     The maximum number of bins
 * @param alpha         The alpha-accuraxy level of q-quantile
 * @return              Parameters of the sketch
 */
extern DDS_type *DDS_Init(int offset, int bin_limit, float alpha);

/**
 * /brief               DDS destructor
 * @param dds           Parameters of the sketch
 */
extern void DDS_Destroy(DDS_type *dds);
/**
 * \brief               Return the number of bins currently in the sketch
 * @param dds           Parameters of the sketch
 * @return              The bins size, (number of bins)
 */
extern long DDS_Size(DDS_type *dds);

/**
 * \brief               Given a value x, getKey returns the bucket index
 * @param dds           Parameters of the sketch
 * @param item          The the input value
 * @return              The index of the bucket containing the value x
 */
extern int DDS_GetKey(DDS_type *dds, double item);

/**
 * \brief               Given a value x, getKey returns the bucket ind
 * @param dds           Parameters of the sketch
 * @param item          The the input value
 * @param ln_gamma      this is not a required parameter; it is defined as log(gamma)
 * @return              The index of the bucket containing the value x
 */
extern int DDS_GetKey(DDS_type *dds, double item, float ln_gamma);

/**
 * \brief               Given a bucket index (i), this function returns the estimation of the rank x_q
 * @param dds           Parameters of the sketch
 * @param i             The key of the bucket
 * @return              The estimate of the rank x_q
 */
extern double DDS_GetRank(DDS_type *dds, int i);

/**
 * \brief               This function creates a new bucket with index associated with the value (item), or if that bucket already exists, it simply add 1 to the bucket's counter
 * @param dds           Parameters of the sketch
 * @param item          The the input value
 * @return              0 success, -1 error
 */
extern int DDS_Add(DDS_type *dds, double item);

/**
 * @brief               This function returns the bound associated with the key (i), (gamma^i)
 * @param dds           Parameters of the sketch
 * @param i             The key of the bucket
 * @return
 */
double DDS_GetBound(DDS_type *dds, int i);

/**
 * @brief               This function returns the bound associated with the key (i), (gamma^i)
 * @param dds           Parameters of the sketch
 * @param i             The key of the bucket
 * @param gamma         Gamma
 * @return
 */
double DDS_GetBound(DDS_type *dds, int i, float gamma);

/**
 * \brief               This function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1 otherwise it simply decrements by 1 the bucket's counter
 * @param dds           Parameters of the sketch
 * @param item          The the input value
 * @return              0 success
 */
extern int DDS_Delete(DDS_type *dds, double item);

/**
 * @brief               In order to reduce the bucket's number, we need to increase the range of the bucket's index.
 * @param dds           Parameters of the sketch
 * @return              0 success, -1 error
 */
extern int DDS_expand(DDS_type *dds);

/**
 * \brief               The function computes the estimate of the desired q-quantile (q)
 * @param dds           Parameters of the sketch
 * @param q             The desired q-quantile
 * @return              The estimate of the desired q-quantile
 */
extern double DDS_GetQuantile(DDS_type *dds, float q);

/**
 * \brief               Merge function merges the bins in dds1 with the bins of dds2 dds1 is the result of the merge operation
 * @param dds1          Parameters of the sketch
 * @param dds2          Parameters of the sketch
 */
extern void DDS_merge(DDS_type *dds1, DDS_type *dds2);

/**
 * @brief               This function computes the sum of all counter of the bins
 * @param dds           Parameters of the sketch
 * @return              int
 */
extern int DDS_SumBins(DDS_type *dds);

/**
 * @brief               This function prints the bins map in a CSV file
 * @param name          File name
 * @param bins          Bins map
 * @return              0 success
 */
extern int DDS_PrintCSV(string name, map<int,int> *bins);

/**
 * @brief               This function checks if all the elements in the stream have a corresponding bucket
 * @param dds           Parameters of the sketch
 * @param item          Input value
 * @return              0 success, -1 failed
 */
extern int DDS_CheckAll(DDS_type *dds, double item);

/**
 * @brief               This function expands all the bins in the map, increasing alpha by 0.01. The values in the old range are redistributed in the new range with a proportional way (supposing all values uniforming distributed on the interval)
 * @param dds           Parameters of the sketch
 * @return              0 success
 */
extern int DDS_expandProportional(DDS_type *dds);

/**
 * @brief               This function collapses the last two buckets
 * @param dds           Parameters of the sketch
 * @return
 */
extern int DDS_Collapse(DDS_type *dds);

/**
 * @brief               This function collapses the adjacent buckets and remap the key in a new hash-map
 * @param dds           Parameters of the sketch
 * @return
 */
extern int DDS_CollapseNeighbors(DDS_type *dds);

/**
 * @brief               This function add an element to the sketch (whit the remapped method)
 * @param dds           Parameters of the sketch
 * @param item          Input value
 * @return
 */
extern int DDS_AddRemapped(DDS_type *dds, double item);

/**
 * @brief               In order to reduce the bucket's number, we need to increase the range of the bucket's index. (with remapped method)
 * @param dds           Parameters of the sketch
 * @param item          Input value
 * @return
 */
extern int DDS_DeleteCollapseNeighborn(DDS_type *dds, double item);

/**
 * @brief               The function collapses the old buckets in the new buckets based on the new range (range ^ 2)
 * @param dds           Parameters of the sketch
 * @return
 */
extern int DDS_CollapsePlus(DDS_type *dds);