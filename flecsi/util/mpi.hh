/*
    @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
   /@@/////  /@@          @@////@@ @@////// /@@
   /@@       /@@  @@@@@  @@    // /@@       /@@
   /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
   /@@////   /@@/@@@@@@@/@@       ////////@@/@@
   /@@       /@@/@@//// //@@    @@       /@@/@@
   /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
   //       ///  //////   //////  ////////  //

   Copyright (c) 2016, Triad National Security, LLC
   All rights reserved.
                                                                              */
#pragma once

/*! @file */

#include <flecsi-config.h>

#include "flecsi/util/serialize.hh"

#if !defined(FLECSI_ENABLE_MPI)
#error FLECSI_ENABLE_MPI not defined! This file depends on MPI!
#endif

#include <complex>
#include <cstddef> // byte
#include <cstdint>
#include <type_traits>

#include <mpi.h>

namespace flecsi {
namespace util {

/*!
 Wrapper to convert from C++ types to MPI types.

 @tparam TYPE The C++ P.O.D. type, e.g., double.

 @ingroup coloring
 */

template<typename TYPE>
struct mpi_typetraits {
  // NB: OpenMPI's predefined handles are not constant expressions.
  template<bool Make = true>
  static MPI_Datatype type() {
    using namespace std;
    static_assert(is_same_v<TYPE, remove_cv_t<remove_reference_t<TYPE>>>);
    // List is from MPI 3.2 draft.
    // TIP: Specializations would collide for, say, int32_t==int.
    // The constexpr-if tree saves on template instantiations.
    if constexpr(is_arithmetic_v<TYPE>) {
      if constexpr(is_same_v<TYPE, char>)
        return MPI_CHAR;
      else if constexpr(is_same_v<TYPE, short>)
        return MPI_SHORT;
      else if constexpr(is_same_v<TYPE, int>)
        return MPI_INT;
      else if constexpr(is_same_v<TYPE, long>)
        return MPI_LONG;
      else if constexpr(is_same_v<TYPE, long long>)
        return MPI_LONG_LONG;
      else if constexpr(is_same_v<TYPE, signed char>)
        return MPI_SIGNED_CHAR;
      else if constexpr(is_same_v<TYPE, unsigned char>)
        return MPI_UNSIGNED_CHAR;
      else if constexpr(is_same_v<TYPE, unsigned short>)
        return MPI_UNSIGNED_SHORT;
      else if constexpr(is_same_v<TYPE, unsigned>)
        return MPI_UNSIGNED;
      else if constexpr(is_same_v<TYPE, unsigned long>)
        return MPI_UNSIGNED_LONG;
      else if constexpr(is_same_v<TYPE, unsigned long long>)
        return MPI_UNSIGNED_LONG_LONG;
      else if constexpr(is_same_v<TYPE, float>)
        return MPI_FLOAT;
      else if constexpr(is_same_v<TYPE, double>)
        return MPI_DOUBLE;
      else if constexpr(is_same_v<TYPE, long double>)
        return MPI_LONG_DOUBLE;
      else if constexpr(is_same_v<TYPE, wchar_t>)
        return MPI_WCHAR;
#ifdef INT8_MIN
      else if constexpr(is_same_v<TYPE, int8_t>)
        return MPI_INT8_T;
      else if constexpr(is_same_v<TYPE, uint8_t>)
        return MPI_UINT8_T;
#endif
#ifdef INT16_MIN
      else if constexpr(is_same_v<TYPE, int16_t>)
        return MPI_INT16_T;
      else if constexpr(is_same_v<TYPE, uint16_t>)
        return MPI_UINT16_T;
#endif
#ifdef INT32_MIN
      else if constexpr(is_same_v<TYPE, int32_t>)
        return MPI_INT32_T;
      else if constexpr(is_same_v<TYPE, uint32_t>)
        return MPI_UINT32_T;
#endif
#ifdef INT64_MIN
      else if constexpr(is_same_v<TYPE, int64_t>)
        return MPI_INT64_T;
      else if constexpr(is_same_v<TYPE, uint64_t>)
        return MPI_UINT64_T;
#endif
      else if constexpr(is_same_v<TYPE, MPI_Aint>)
        return MPI_AINT;
      else if constexpr(is_same_v<TYPE, MPI_Offset>)
        return MPI_OFFSET;
      else if constexpr(is_same_v<TYPE, MPI_Count>)
        return MPI_COUNT;
      else if constexpr(is_same_v<TYPE, bool>)
        return MPI_CXX_BOOL;
      else
        return make<Make>();
    }
    else if constexpr(is_same_v<TYPE, complex<float>>)
      return MPI_CXX_FLOAT_COMPLEX;
    else if constexpr(is_same_v<TYPE, complex<double>>)
      return MPI_CXX_DOUBLE_COMPLEX;
    else if constexpr(is_same_v<TYPE, complex<long double>>)
      return MPI_CXX_LONG_DOUBLE_COMPLEX;
    else if constexpr(is_same_v<TYPE, byte>)
      return MPI_BYTE;
    else
      return make<Make>();
  }

private:
  template<bool Make>
  static MPI_Datatype make() {
    static_assert(Make, "type not predefined");
    static_assert(std::is_trivially_copyable_v<TYPE>);
    // TODO: destroy at MPI_Finalize
    static const MPI_Datatype ret = [] {
      MPI_Datatype data_type;
      MPI_Type_contiguous(sizeof(TYPE), MPI_BYTE, &data_type);
      MPI_Type_commit(&data_type);
      return data_type;
    }();
    return ret;
  }
};

template<class T>
MPI_Datatype
mpi_type() {
  return mpi_typetraits<T>::type();
}
// Use references to this to restrict to predefined types before MPI_Init:
template<class T>
extern const MPI_Datatype // 'extern' works around GCC bug #90493
  mpi_static_type = mpi_typetraits<T>::template type<false>();

namespace mpi {

// Convenience variables
inline const auto size_type = mpi_typetraits<std::size_t>::type();
inline const auto byte_type = mpi_typetraits<std::byte>::type();
inline const auto char_type = mpi_typetraits<char>::type();
inline const auto unsigned_char_type = mpi_typetraits<unsigned char>::type();
inline const auto short_type = mpi_typetraits<short>::type();
inline const auto unsigned_short_type = mpi_typetraits<unsigned short>::type();
inline const auto int_type = mpi_typetraits<int>::type();
inline const auto unsigend_type = mpi_typetraits<unsigned>::type();
inline const auto long_type = mpi_typetraits<long>::type();
inline const auto float_type = mpi_typetraits<float>::type();
inline const auto double_type = mpi_typetraits<double>::type();
inline const auto long_double_type = mpi_typetraits<long double>::type();

/*!
  Convenience function to get basic MPI communicator information.
 */

inline auto
info(MPI_Comm comm = MPI_COMM_WORLD) {
  int rank, size;
  MPI_Group group;

  MPI_Comm_size(comm, &size);
  MPI_Comm_group(comm, &group);
  MPI_Group_rank(group, &rank);

  return std::make_tuple(rank, size, group);
} // info

/*!
  One-to-All (variable) communication pattern.

  This function uses the FleCSI serialization interface with a packing functor
  to communicate data from the root rank (0) to all other ranks.

  @tparam F The packing functor type, which must define a \em return_type,
            and the \emph () operator, taking \em rank and \em size as integer
            arguments. The \em return_type is the type return by the packing
            functor.

  @param f    An instance of the packing functor.
  @param comm An MPI communicator.

  @return For all ranks besides the root rank (0), the communicated data. For
          the root rank (0), the functor applied for the root rank (0) and size.
 */

template<typename F>
inline auto
one_to_allv(F const & f, MPI_Comm comm = MPI_COMM_WORLD) {
  auto [rank, size, group] = info(comm);

  std::vector<MPI_Request> requests;

  if(rank == 0) {
    for(size_t r{1}; r < std::size_t(size); ++r) {
      std::vector<std::byte> data(serial_put(f(r, size)));
      const std::size_t bytes = data.size();

      requests.resize(requests.size() + 1);
      MPI_Isend(&bytes, 1, mpi::size_type, r, 0, comm, &requests.back());
      requests.resize(requests.size() + 1);
      MPI_Isend(
        data.data(), bytes, mpi::byte_type, r, 0, comm, &requests.back());
    } // for

    std::vector<MPI_Status> status(requests.size());
    MPI_Waitall(requests.size(), requests.data(), status.data());
  }
  else {
    std::size_t bytes{0};
    MPI_Status status;
    MPI_Recv(&bytes, 1, mpi::size_type, 0, 0, comm, &status);
    std::vector<std::byte> data(bytes);
    MPI_Recv(data.data(), bytes, mpi::byte_type, 0, 0, comm, &status);
    auto const * p = data.data();
    return serial_get<typename F::return_type>(p);
  } // if

  return f(0, size);
} // one_to_allv

/*!
  All-to-All (variable) communication pattern implemented with non-blocking
  send and receive operations.

  This function uses the FleCSI serialization interface with a packing functor
  to communicate data from all ranks to all other ranks.

  @tparam F The packing functor type, which must define a \em return_type,
            a \em count() method, taking the calling rank as an argument,
            and the \emph () operator, taking \em rank and \em size as integer
            arguments. The \em return_type is the type return by the packing
            functor.

  @param f    An instance of the packing functor.
  @param comm An MPI communicator.

  @return A std::vector<return_type>.
 */

template<typename F>
inline auto
all_to_allv(F const & f, MPI_Comm comm = MPI_COMM_WORLD) {
  auto [rank, size, group] = info(comm);

  std::vector<std::size_t> counts;
  counts.reserve(size);
  std::vector<std::size_t> bytes(size);

  for(std::size_t r{0}; r < std::size_t(size); ++r) {
    counts.emplace_back(f.count(r));
  } // for

  MPI_Alltoall(
    counts.data(), 1, mpi::size_type, bytes.data(), 1, mpi::size_type, comm);

  std::vector<MPI_Request> requests;

  std::vector<std::vector<std::byte>> send_bufs(size);
  std::vector<std::vector<std::byte>> recv_bufs(size);

  for(std::size_t r{0}; r < std::size_t(size); ++r) {
    if(bytes[r] > 0) {
      recv_bufs[r].resize(bytes[r]);
      requests.resize(requests.size() + 1);
      MPI_Irecv(recv_bufs[r].data(),
        recv_bufs[r].size(),
        mpi::byte_type,
        r,
        0,
        comm,
        &requests.back());
    } // if
  } // for

  for(std::size_t r{0}; r < std::size_t(size); ++r) {
    if(counts[r] > 0) {
      requests.resize(requests.size() + 1);
      send_bufs[r] = serial_put(f(r, size));
      MPI_Isend(send_bufs[r].data(),
        send_bufs[r].size(),
        mpi::byte_type,
        r,
        0,
        comm,
        &requests.back());
    } // if
  } // for

  std::vector<MPI_Status> status(requests.size());
  MPI_Waitall(requests.size(), requests.data(), status.data());

  std::vector<typename F::return_type> result;
  result.reserve(size);
  for(std::size_t r{0}; r < std::size_t(size); ++r) {
    if(recv_bufs[r].size() > 0) {
      auto const * p = recv_bufs[r].data();
      result.emplace_back(serial_get<typename F::return_type>(p));
    }
  } // for

  return result;
} // all_to_allv

} // namespace mpi

} // namespace util
} // namespace flecsi
