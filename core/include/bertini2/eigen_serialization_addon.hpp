


// the following code comes from
// https://stackoverflow.com/questions/18382457/eigen-and-boostserialize
// and adds support for serialization of Eigen types
//
// question asked by user Gabriel and answered by iNFINITEi
//
// answer code repeated here verbatim.
// please update this comment if this code is changed,
// and post the modifications to the above referenced post on SO.



/**
 * @file eigen_serialization_addon.hpp
 */

#ifndef EIGEN_SERIALIZATION_ADDON_HPP
#define EIGEN_SERIALIZATION_ADDON_HPP

#pragma once




friend class boost::serialization::access;
template<class Archive>
void save(Archive & ar, const unsigned int version) const {
  derived().eval();
  const Eigen::Index rows = derived().rows(), cols = derived().cols();
  ar & rows;
  ar & cols;
  for (Index j = 0; j < cols; ++j )
    for (Index i = 0; i < rows; ++i )
      ar & derived().coeff(i, j);
}

template<class Archive>
void load(Archive & ar, const unsigned int version) {
  Eigen::Index rows, cols;
  ar & rows;
  ar & cols;
  if (rows != derived().rows() || cols != derived().cols() )
    derived().resize(rows, cols);
  ar & boost::serialization::make_array(derived().data(), derived().size());
}

template<class Archive>
void serialize(Archive & ar, const unsigned int file_version) {
  boost::serialization::split_member(ar, *this, file_version);
}

#endif // EIGEN_SERIALIZATION_ADDON_HPP
