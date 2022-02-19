#include "header.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


namespace py = pybind11;


PYBIND11_MODULE(abp, m) {
    m.doc() = "pybind11 of Active Brownian Particle dynamics"; 

    // Class point
    py::class_<point>(m, "point")
        .def_readwrite("x", &point::x)
        .def_readwrite("y", &point::y)
        .def(py::init<double, double>(), py::arg("x"), py::arg("y"))
        .def("distance_to_point", &point::distance_to_point,  py::arg("position"));
        
    // Derived class region
    py::class_<region, point>(m, "region")
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"),  py::arg("radius"))
        .def_readwrite("radius", &region::radius);

    m.def("getVector", &getVector, "Compute the vector(point) connecting two points",  py::arg("end"),  py::arg("start"));

    // Class active brownian particle
    py::class_<ABP_2d>(m, "ABP_2d")
        // Static members
        .def_readwrite("position_x", &ABP_2d::position_x)
        .def_readwrite("position_y", &ABP_2d::position_y)
        .def_readwrite("theta", &ABP_2d::theta)
        .def_readwrite("bool_reactant", &ABP_2d::bool_reactant)
        .def_readwrite("bool_target", &ABP_2d::bool_target)
        .def_readwrite("reactive_path", &ABP_2d::reactive_path)
        .def_readwrite("num_steps", &ABP_2d::num_steps)
        .def_readwrite("dt", &ABP_2d::dt)
        .def_readwrite("v", &ABP_2d::v)
        .def_readwrite("D_r", &ABP_2d::D_r)
        .def_readwrite("D_theta", &ABP_2d::D_theta)
        .def_readwrite("k", &ABP_2d::k)
        .def_readwrite("L", &ABP_2d::L)
        .def_readwrite("mu", &ABP_2d::mu)
        .def_readwrite("w", &ABP_2d::w)
        .def_readwrite("reactant", &ABP_2d::reactant)
        // Functions
        .def(py::init<const region&,const region&, unsigned&, double&, double&, double&, double&, double&, double&, double&, double&>(),
        py::arg("reactant"), py::arg("target"), py::arg("num_steps"), py::arg("dt"), py::arg("v"), py::arg("D_r"), py::arg("D_theta"),
        py::arg("k"), py::arg("L"), py::arg("mu"), py::arg("w"))
        .def("apply_pbc", &ABP_2d::apply_pbc, py::arg("position"))
        .def("pbc_distance", &ABP_2d::pbc_distance, py::arg("A"), py::arg("B"))
        .def("potential", &ABP_2d::potential, py::arg("position"))
        .def("compute_force", &ABP_2d::compute_force, py::arg("position"))
        .def("position_step", &ABP_2d::position_step, py::arg("position"), py::arg("theta"), py::arg("noise_x"), py::arg("noise_y"))
        .def("theta_step", &ABP_2d::theta_step, py::arg("theta"), py::arg("noise_theta"))
        .def("is_near_minimum", &ABP_2d::is_near_minimum, py::arg("position"))
        .def("is_inside_region", &ABP_2d::is_inside_region, py::arg("position"), py::arg("target"))
        .def("dynamics", &ABP_2d::dynamics)
        .def("print_dynamics", &ABP_2d::print_dynamics, py::arg("filename"))
        .def("print_bool_dynamics", &ABP_2d::print_bool_dynamics, py::arg("filename"));
}
