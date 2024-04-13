#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "_wrapper.h"

namespace py = pybind11;

using c_array = py::array_t<double, py::array::c_style | py::array::forcecast>;

struct TM_result {
  c_array t; // optimal translation
  c_array u; // optimal rotation

  double TM1, TM2; // normalized TM scores

  TM_result(const double t_[3], const double u_[3][3], double TM1_, double TM2_)
      : t(std::vector<ptrdiff_t>{3}), u(std::vector<ptrdiff_t>{3, 3}),
        TM1(TM1_), TM2(TM2_) {
    auto r_t = t.mutable_unchecked<1>();
    auto r_u = u.mutable_unchecked<2>();

    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        r_u(i, j) = u_[i][j];
      }
      r_t[i] = t_[i];
    }
  }
};

static void _check_shape(const c_array &data, std::string name, size_t which,
                         size_t expected) {
  if (data.shape(which) != expected) {
    std::stringstream s;
    s << "Incorrect shape " << which << " for array '" << name << "' "
      << "(expected " << expected << " but got " << data.shape(which) << ")";
    throw std::runtime_error(s.str());
  }
}

static std::vector<double *> _to_raw(const c_array &arr) {
  py::buffer_info buf = arr.request();
  double *ptr = static_cast<double *>(buf.ptr);

  std::vector<double *> data(arr.shape(0));

  for (size_t i = 0; i < arr.shape(0); i++) {
    data[i] = ptr;
    ptr += arr.shape(1);
  }

  return data;
}

static TM_result tm_align(c_array x, c_array y, std::string seqx,
                          std::string seqy) {
  // coordinates
  _check_shape(x, "x", 0, seqx.size());
  _check_shape(y, "y", 0, seqy.size());
  _check_shape(x, "x", 1, 3);
  _check_shape(y, "y", 1, 3);

  auto raw_x = _to_raw(x);
  auto raw_y = _to_raw(y);

  // output parameters
  double TM1, TM2;
  double t[3];
  double u[3][3];

  _tmalign_wrapper(raw_x.data(), raw_y.data(), seqx.c_str(), seqy.c_str(),
                   seqx.size(), seqy.size(), t, u, TM1, TM2);

  return TM_result(t, u, TM1, TM2);
}

const char* tm_align_docstring =
    "Perform structural aligment of proteins using the TM-score algorithm.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "x : ndarray\n"
    "    Coordinates of the backbone of protein 1\n"
    "y : ndarray\n"
    "    Coordinates of the backbone of protein 2\n"
    "seqx : str\n"
    "    Protein sequence 1\n"
    "seqy : str\n"
    "    Protein sequence 2\n"
    "\n"
    "Returns\n"
    "-------\n"
    "result : TM_result\n"
    "    An object with the optimal translation and rotation, as well as\n"
    "    the TM-scores.\n";


PYBIND11_MODULE(_bindings, m) {
  m.doc() = "Low-level Python wrappers for tm-align";
  m.def("tm_align",
        &tm_align,
        tm_align_docstring,
        py::arg("x"),
        py::arg("y"),
        py::arg("seqx"),
        py::arg("seqy"));

  py::class_<TM_result>(m, "TMResult", "Results wrapper for the TM-align algorithm")
      .def_readonly("t", &TM_result::t, "Optimal translation from protein 1 to 2")
      .def_readonly("u", &TM_result::u, "Optimal rotation from protein 1 to 2")
      .def_readonly("tm_norm_chain2", &TM_result::TM1,
                    "TM-score normalized by the length of protein 2")
      .def_readonly("tm_norm_chain1", &TM_result::TM2,
                    "TM-score normalized by the length of protein 1");
}
