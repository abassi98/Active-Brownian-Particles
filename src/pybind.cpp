#include "header.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


namespace py = pybind11;


PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 of Active Brownian Particle dynamics"; 

    // Class point
    py::class_<point>(m, "point")
        .def_readwrite("x", &point::x)
        .def_readwrite("y", &point::y)
        .def(py::init<double, double>())
        .def("distance_to_point", &point::distance_to_point);
        
    // Derived class region
    py::class_<region, point>(m, "region")
        .def(py::init<double, double, double>())
        .def_readwrite("radius", &region::radius);

    m.def("getVector", &getVector, "Compute the vector(point) connecting two points");

    // Class active brownian particle
    py::class_<ABP_2d>(m, "ABP_2d")
        // Static members
        .def_readwrite("position_x", &ABP_2d::position_x)
        .def_readwrite("position_y", &ABP_2d::position_y)
        .def_readwrite("theta", &ABP_2d::theta)
        .def_readwrite("bool_reactant", &ABP_2d::bool_reactant)
        .def_readwrite("bool_target", &ABP_2d::bool_target)
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
        .def(py::init<const region&,const region&, unsigned&, double&, double&, double&, double&, double&, double&, double&, double&>())
        .def("apply_pbc", &ABP_2d::apply_pbc)
        .def("pbc_distance", &ABP_2d::pbc_distance)
        .def("potential", &ABP_2d::potential)
        .def("compute_force", &ABP_2d::compute_force)
        .def("position_step", &ABP_2d::position_step)
        .def("theta_step", &ABP_2d::theta_step)
        .def("is_near_minimum", &ABP_2d::is_near_minimum)
        .def("is_inside_region", &ABP_2d::is_inside_region)
        .def("dynamics", &ABP_2d::dynamics)
        .def("print_dynamics", &ABP_2d::print_dynamics)
        .def("print_bool_dynamics", &ABP_2d::print_bool_dynamics);
}
