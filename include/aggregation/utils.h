//
// Created by Kazem on 10/11/19.
//

#ifndef PROJECT_UTILS_H
#define PROJECT_UTILS_H

#include <vector>
#include "def.h"

namespace sym_lib{

/// Comparing two time struct
/// \param a
/// \param b
/// \return
 bool time_cmp(timing_measurement a, timing_measurement b);


 /// Finding the median time
 /// \param time_array
 /// \return
 timing_measurement
 time_median(std::vector<timing_measurement> time_array);


/// Safely compute a*k, where k should be small, and check for integer overflow.
/// If overflow occurs, return 0 and set OK to FALSE.  Also return 0 if OK is
/// FALSE on input.
/// \param a
/// \param k
/// \param ok
/// \return
size_t mult_size_t(size_t a, size_t k, int *ok);

/// Safely compute a+b, and check for integer overflow.
/// If overflow occurs, return 0 and set OK to FALSE.
/// Also return 0 if OK is FALSE on input.
/// \param a
/// \param b
/// \param ok
/// \return
 size_t add_size_t (size_t a, size_t b, int *ok);


/// Partitions a flat set based on the given weight
/// \param n
/// \param set : input set
/// \param weight : is the weight[i] of node in set[i]
/// \param n_parts : number of partitions requested
/// \param target_weight : the weight of each requested partition.
/// \return
 int *partition_by_weight(int n,const int *set, const double *weight,
                          int n_parts,
                          double *target_weight);


 /// partitions a DAG by doing bfs till reachs the weight threshold
 /// \param n in
 /// \param df in
 /// \param weight in
 /// \param final_level_no  out
 /// \param fina_level_ptr out
 /// \param part_no in
 /// \param final_part_ptr  out
 /// \param final_node_ptr out
 /// \param target_weight in
 void partition_by_bfs(int n, CSC *df, const double *weight,
                       int &final_level_no, // will be one in this case
                       int* &fina_level_ptr, int part_no,
                       int *&final_part_ptr, int *&final_node_ptr,
                       double *target_weight);


 /// Calculates the summation of a vector
 /// \tparam type
 /// \param n
 /// \param vec
 /// \return
 template<class type> type sum_vector(int n, type *vec){
  double sum = 0;
  for (int i = 0; i < n; ++i) sum += vec[i];
  return sum;
 }



 // trim from start (in place)
 static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
   return !std::isspace(ch);
  }));
 }

// trim from end (in place)
 static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
   return !std::isspace(ch);
  }).base(), s.end());
 }

// trim from both ends (in place)
 static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
 }


 /// find the name of a file from a path
 std::string strip_name(std::string name);

}
#endif //PROJECT_UTILS_H
