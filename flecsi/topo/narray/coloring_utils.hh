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

#if !defined(__FLECSI_PRIVATE__)
#error Do not include this file directly!
#endif

#include "flecsi/flog.hh"
#include "flecsi/topo/narray/types.hh"
#include "flecsi/util/color_map.hh"

#include <functional>
#include <optional>
#include <vector>

namespace flecsi {
namespace topo {
namespace narray_impl {

inline std::vector<std::size_t>
sieve(std::size_t n) {
  std::vector<bool> prime(n + 1, true);
  std::vector<std::size_t> ret;

  for(std::size_t p{2}; p * p <= n; ++p) {
    if(prime[p]) {
      for(std::size_t i{2 * p}; i <= n; i += p) {
        prime[i] = false;
      }
    }
  }

  for(std::size_t p{2}; p <= n; ++p) {
    if(prime[p])
      ret.push_back(p);
  }

  return ret;
} // sieve

inline std::vector<std::size_t>
factor(std::size_t np) {
  std::vector<std::size_t> facs;
  auto primes = sieve(np);

  std::size_t p{0};
  while(np != 1) {
    if(np % primes[p] == 0) {
      facs.push_back(primes[p]);
      np = np / primes[p];
    }
    else {
      ++p;
    }
  }

  std::sort(facs.begin(), facs.end());
  std::reverse(facs.begin(), facs.end());

  return facs;
} // factor

std::vector<std::size_t>
distribute(std::size_t np, std::vector<std::size_t> indices) {
  std::vector<std::size_t> parts(indices.size(), 1);

  auto facs = factor(np);

  // greedy decomposition
  for(auto fac : facs) {
    auto maxind = std::distance(
      indices.begin(), std::max_element(indices.begin(), indices.end()));
    parts[maxind] *= fac;
    indices[maxind] /= fac;
  }

  return parts;
} // decomp

/*
  Using the index colors for the given color, this function determines the
  position of the local partition within the color space. The "faces"
  variable is a bit-array that stores (for each axis) boolean values
  for "low", and "high". If an axis is neither low nor high, it is
  interior. The enumeration defining these masks is in types.hh.
 */

auto
orientation(std::size_t dimension,
  std::vector<std::size_t> const & color_indices,
  std::vector<std::size_t> const & axis_colors) {
#define FACE(e, c, c0, c1, c2)                                                 \
  ((((e) == 0) && ((e) == ((c)-1)))                                            \
      ? c1 | c2                                                                \
      : ((0 < (e) && e < ((c)-1)) ? c0 : ((e) == 0 ? c1 : c2)))

  std::uint32_t faces = interior;
  std::uint32_t shft{low};
  for(std::size_t axis{0}; axis < dimension; ++axis) {
    faces |=
      FACE(color_indices[axis], axis_colors[axis], interior, shft, shft << 1);
    shft <<= 2;
  } // for
#undef FACE

  return faces;
} // orientation

/*
 */

template<std::size_t HaloDepth, std::size_t BoundaryDepth>
inline auto
make_color(std::size_t dimension,
  std::vector<std::size_t> const & color_indices,
  std::vector<util::color_map> const & axcm,
  uint32_t faces) {
  index_coloring idxco;
  idxco.extents.resize(dimension);
  idxco.global.resize(dimension);
  idxco.owned[0].resize(dimension);
  idxco.owned[1].resize(dimension);
  idxco.exclusive[0] = idxco.owned[0];
  idxco.exclusive[1] = idxco.owned[1];

  std::vector</* over axes */
    std::vector</* over intervals */
      std::pair<std::size_t, /* owner color */
        std::pair<std::size_t, std::size_t> /* interval */
        >>>
    ghstitvls(dimension);

  for(std::size_t axis{0}; axis < dimension; ++axis) {
    auto axis_color = color_indices[axis];
    idxco.global[axis] = axcm[axis].index_offset(color_indices[axis], 0);

    idxco.extents[axis] = axcm[axis].indices(axis_color, 0);
    idxco.owned[0][axis] = 0;
    idxco.owned[1][axis] = axcm[axis].indices(axis_color, 0);

    std::uint32_t bits = faces >> axis * 2;

    if(bits & low && bits & high) {
      /*
        This is a degenerate dimension, i.e., it is flat with a single
        color layer. Therefore, we do not add halo extensions.
       */

      idxco.extents[axis] += 2 * BoundaryDepth;

      idxco.owned[0][axis] += BoundaryDepth;
      idxco.owned[1][axis] += BoundaryDepth;

      idxco.exclusive[0][axis] = idxco.owned[0][axis];
      idxco.exclusive[1][axis] = idxco.owned[1][axis];
    }
    else if(bits & low) {
      /*
        This dimension is a low edge.
       */

      idxco.extents[axis] += BoundaryDepth + HaloDepth;

      idxco.owned[0][axis] += BoundaryDepth;
      idxco.owned[1][axis] += BoundaryDepth;

      idxco.exclusive[0][axis] = idxco.owned[0][axis];
      idxco.exclusive[1][axis] = idxco.owned[1][axis] - HaloDepth;

      ghstitvls[axis].push_back({color_indices[axis] + 1,
        {idxco.owned[1][axis], idxco.owned[1][axis] + HaloDepth}});
    }
    else if(bits & high) {
      /*
        This dimension is a high edge.
       */

      idxco.extents[axis] += HaloDepth + BoundaryDepth;

      idxco.global[axis] -= HaloDepth;

      idxco.owned[0][axis] += HaloDepth;
      idxco.owned[1][axis] += HaloDepth;

      idxco.exclusive[0][axis] = idxco.owned[0][axis] + HaloDepth;
      idxco.exclusive[1][axis] = idxco.owned[1][axis];

      ghstitvls[axis].push_back({color_indices[axis] - 1,
        {idxco.owned[0][axis] - HaloDepth, idxco.owned[0][axis]}});
    }
    else {
      /*
        This dimension is interior.
       */

      idxco.extents[axis] += 2 * HaloDepth;

      idxco.global[axis] -= HaloDepth;

      idxco.owned[0][axis] += HaloDepth;
      idxco.owned[1][axis] += HaloDepth;

      idxco.exclusive[0][axis] = idxco.owned[0][axis] + HaloDepth;
      idxco.exclusive[1][axis] = idxco.owned[1][axis] - HaloDepth;

      ghstitvls[axis].push_back({color_indices[axis] + 1,
        {idxco.owned[1][axis], idxco.owned[1][axis] + HaloDepth}});
      ghstitvls[axis].push_back({color_indices[axis] - 1,
        {idxco.owned[0][axis] - HaloDepth, idxco.owned[0][axis]}});
    } // if
  } // for

  return std::make_pair(idxco, ghstitvls);
} // make_color

/*
 */

struct coloring_definition {
  coord axis_colors;
  coord axis_extents;
  bool diagonals = false;
}; // struct coloring_definition

template<std::size_t HaloDepth, std::size_t BoundaryDepth>
inline auto
color(std::vector<coloring_definition> const & index_spaces,
  MPI_Comm comm = MPI_COMM_WORLD) {

  flog_assert(index_spaces.size() != 0, "no index spaces defined");

  const std::size_t dimension = index_spaces[0].axis_colors.size();

  flog_assert(dimension < 17,
    "current implementation is limited to 16 dimensions (uint32_t)");

  std::size_t colors{1};
  auto [rank, size] = util::mpi::info(comm);
  std::vector<std::map<std::size_t, index_coloring>> colorings;

  for(std::size_t is{0}; is < index_spaces.size(); ++is) {
    auto const axis_colors = index_spaces[is].axis_colors;
    auto const axis_extents = index_spaces[is].axis_extents;
    bool const diagonals = index_spaces[is].diagonals;

    flog_assert(axis_colors.size() == axis_extents.size(),
      "argument mismatch: sizes(" << axis_colors.size() << "vs. "
                                  << axis_extents.size()
                                  << ") must be consistent");
    flog_assert(axis_colors.size() == dimension,
      "size must match the intended dimension(" << dimension << ")");

    /*
      Create a color map for each dimension. Because we are using a
      tensor-product strategy below, we can use each map to define the
      sub colors and offsets for an individual axis.
     */

    std::size_t idx_colors{1};
    std::size_t indices{1};
    std::vector<util::color_map> axcm;
    for(std::size_t d{0}; d < dimension; ++d) {
      if(is == 0) {
        colors *= axis_colors[d];
      }
      idx_colors *= axis_colors[d];
      indices *= axis_extents[d];
      axcm.emplace_back(
        util::color_map{axis_colors[d], axis_colors[d], axis_extents[d]});
    } // for

    flog_assert(idx_colors == colors,
      "index spaces must have a consistent number of colors");

    /*
      Create a color map for the total number of colors (product of axis
      colors) to the number of processes.
     */

    util::color_map cm(size, idx_colors, indices);

    /*
      Create a coloring for each color on this process.
     */

    flog(warn) << "rank colors: " << cm.colors(rank) << std::endl;
    std::map<std::size_t, index_coloring> coloring;
    for(std::size_t c{0}; c < cm.colors(rank); ++c) {
      /*
        Convenience functions to map between colors and indices.
       */

      auto co2idx = [](std::size_t co, coord szs) {
        coord indices;
        for(auto sz : szs) {
          indices.emplace_back(co % sz);
          co /= sz;
        }
        return indices;
      };

      auto idx2co = [](coord const & idx, coord const & szs) {
        std::size_t co{0}, pr{szs[0]};
        for(std::size_t i{0}; i < idx.size() - 1; ++i) {
          co += idx[i + 1] * pr;
          pr *= szs[i + 1];
        }
        return co + idx[0];
      };

      auto color_indices = co2idx(cm.color_id(rank, c), axis_colors);
#if 0
      flog(warn) << "color: " << cm.color_offset(rank) + c << std::endl;
      flog(warn) << log::insert(color_indices, "indices") << std::endl;
      flog(warn) << "color check: " << idx2co(color_indices, axis_colors)
                 << std::endl;
#endif

      /*
        Find our orientation within the color space.
       */

      uint32_t faces = orientation(dimension, color_indices, axis_colors);

      /*
        Make the coloring information for our color.
       */

      auto [idxco, ghstitvls] = make_color<HaloDepth, BoundaryDepth>(
        dimension, color_indices, axcm, faces);

      // This seems like a compiler bug!
      auto idx_coloring = idxco;
      auto ghost_intervals = ghstitvls;

      /*
        This is the "subsequent" loop referenced above. Here we compose
        the intervals from each sub-dimension to form the actual
        full-dimensional subregions. These define the coloring.
       */

      std::function<std::vector<std::pair<coord, rect>>(std::size_t)> expand =
        [idx_coloring, color_indices, ghost_intervals, diagonals, &expand](
          std::size_t dim) {
          std::vector<std::pair<coord, rect>> sregs;

          for(std::size_t axis{0}; axis < dim; ++axis) {
            if(sregs.size()) {
              /*
                Expand the subregions from the lower dimensions.
               */

              auto subs = sregs;
              sregs.clear();
              for(size_t off{idx_coloring.owned[0][axis]};
                  off < idx_coloring.owned[1][axis];
                  ++off) {
                for(auto s : subs) {
                  s.first[axis] = color_indices[axis];
                  s.second[0][axis] = off;
                  s.second[1][axis] = off + 1;
                  sregs.emplace_back(s);
                } // for
              } // for
            } // if

            /*
              Add the subregions for this dimension.
             */

            for(auto i : ghost_intervals[axis]) {
              coord co(dim, 0);
              coord start(dim, 0);
              coord end(dim, 0);
              for(std::size_t a{0}; a < axis; ++a) {
                co[a] = color_indices[a];
                start[a] = idx_coloring.owned[0][a];
                end[a] = idx_coloring.owned[1][a];
              } // for

              co[axis] = i.first;
              start[axis] = i.second.first;
              end[axis] = i.second.second;
              sregs.push_back({co, rect{start, end}});

              flog_assert(!diagonals, "diagonal support is not implemented");
              if(diagonals) {
                /*
                  Recurse to pull up lower dimensions.
                 */

                auto ssubs = expand(dim - 1);

                /*
                  Add axis information from this dimension to new diagonals.
                 */

                for(auto ss : ssubs) {
                  ss.first[axis] = i.first;
                  ss.second[0][axis] = i.second.first;
                  ss.second[1][axis] = i.second.second;
                  sregs.emplace_back(ss);
                } // for
              } // if
            } // for
          } // for

          return sregs;
        };

      auto subregions = expand(dimension);

#if 0
      for(auto s : subregions) {
        std::stringstream ss;
        ss << "start: [";
        for(std::size_t axis{0}; axis < dimension; ++axis) {
          ss << s.second[0][axis];
          if(axis != dimension - 1)
            ss << ", ";
        } // for
        ss << "] end: (";
        for(std::size_t axis{0}; axis < dimension; ++axis) {
          ss << s.second[1][axis];
          if(axis != dimension - 1)
            ss << ", ";
        } // for
        ss << ") color: " << idx2co(s.first, axis_colors);
        flog(warn) << ss.str() << std::endl;
      } // for
#endif

      /*
        Compute a remote index from a global coordinate.
       */

      auto rmtidx = [dimension](
                      index_coloring const & idxco, coord const & gidx) {
        coord result(dimension);
        for(std::size_t axis{0}; axis < dimension; ++axis) {
          result[axis] = gidx[axis] - idxco.global[axis];
        }
        return result;
      };

      /*
        Map a local coordinate to a global one.
       */

      auto l2g = [dimension](index_coloring const & idxco, coord const & idx) {
        coord result(dimension);
        for(std::size_t axis{0}; axis < dimension; ++axis) {
          result[axis] = idxco.global[axis] + idx[axis];
        }
        return result;
      };

      /*
        The intervals computed in the tensor product strategy above are
        closed on the start of the interval, and open on the end. This
        function is used below to close the end, so that the interval
        can be converted into a memory offset interval.
       */

      auto op2cls = [dimension](coord const & idx) {
        coord result(dimension);
        for(std::size_t axis{0}; axis < dimension; ++axis) {
          result[axis] = idx[axis] - 1;
        }
        return result;
      };

      /*
        Loop through the subregions and create the actual coloring.
       */

      std::unordered_map<std::size_t, index_coloring> idxmap;
      for(auto s : subregions) {
        auto co = idx2co(s.first, axis_colors);
        if(idxmap.find(co) == idxmap.end()) {
          /*
            Create basic coloring information for the owning color, so
            that we can determine the remote offsets for our points.
            The coloring informaiton is stored for subsequent use.
           */

          auto [ridxco, rghstitvls] =
            make_color<HaloDepth, BoundaryDepth>(dimension,
              s.first,
              axcm,
              orientation(dimension, s.first, axis_colors));
          idxmap[co] = ridxco;
        } // if

        // Compute the local memory interval.
        auto const end = idx2co(op2cls(s.second[1]), idxco.extents);
        auto const start = idx2co(s.second[0], idxco.extents);

        // The output intervals are closed on the start
        // and open on the end, i.e., [start, end)
        idxco.intervals.push_back({start, end + 1});

        /*
          Loop through the local interval sizes, and add the remote pointer
          offsets.
         */

        auto const gidx = l2g(idxco, s.second[0]);
        auto const ridx = rmtidx(idxmap.at(co), gidx);
        auto rmtoff = idx2co(ridx, idxmap.at(co).extents);

        for(std::size_t off{0}; off < (end + 1) - start; ++off) {
          idxco.points[co].push_back({start + off, rmtoff + off});
        } // for
      } // for

#if 1
      for(auto i : idxco.intervals) {
        flog(warn) << "<" << i.first << "," << i.second << ">" << std::endl;
      } // for

      for(auto i : idxco.points) {
        std::stringstream ss;
        ss << "color: " << i.first << std::endl;
        for(auto e : i.second) {
          ss << "<" << e.first << "," << e.second << ">" << std::endl;
        } // for
        flog(warn) << ss.str() << std::endl;
      } // for
#endif
      coloring.emplace(c, idxco);
    } // for

    colorings.emplace_back(coloring);
  } // for

  return std::make_pair(colors, colorings);
} // color

} // namespace narray_impl
} // namespace topo
} // namespace flecsi
