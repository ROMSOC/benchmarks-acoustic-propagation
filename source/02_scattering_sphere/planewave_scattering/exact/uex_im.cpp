/* ------------------------------------------------------------------*\
                        ╦═╗╔═╗╔╦╗╔═╗╔═╗╔═╗
                        ╠╦╝║ ║║║║╚═╗║ ║║  
                        ╩╚═╚═╝╩ ╩╚═╝╚═╝╚═╝
 Reduced Order Modelling, Simulation, Optimization of Coupled Systems 
                            2017-2021

 Authors : 
 Ashwin Nayak, Andres Prieto, Daniel Fernandez Comesana
 
 Disclaimer :
 In downloading this SOFTWARE you are deemed to have read and agreed to 
 the following terms: This SOFTWARE has been designed with an exclusive 
 focus on civil applications. It is not to be used for any illegal, 
 deceptive, misleading or unethical purpose or in any military 
 applications. This includes ANY APPLICATION WHERE THE USE OF THE 
 SOFTWARE MAY RESULT IN DEATH, PERSONAL INJURY OR SEVERE PHYSICAL 
 OR ENVIRONMENTAL DAMAGE. Any redistribution of the software must 
 retain this disclaimer. BY INSTALLING, COPYING, OR OTHERWISE USING 
 THE SOFTWARE, YOU AGREE TO THE TERMS ABOVE. IF YOU DO NOT AGREE TO 
 THESE TERMS, DO NOT INSTALL OR USE THE SOFTWARE

 Acknowledgements:
 The ROMSOC project has received funding from the European Union’s 
 Horizon 2020 research and innovation programme under the Marie 
 Skłodowska-Curie Grant Agreement No. 765374.
\*-------------------------------------------------------------------*/

// Include necessary PYBIND files
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

// Include necessary DOLFIN classes
#include <dolfin/function/Expression.h>
#include <dolfin/function/Constant.h>

// Include necessary BOOST classes
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/hankel.hpp>

#include <complex>
#include <iostream>
class ScatteringExact_Im : public dolfin::Expression {
public:
double p0, k, a;
ScatteringExact_Im() : dolfin::Expression(), p0(1.0), k(400.),  a(0.05){}

void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const {

    using namespace boost::math;

    std::complex<double> one_i(0.0,1.0);

    double r = x.norm();

    // Compute Scattered field
    double tol = 1e-10;
    int l=0;
    int N=50;
    std::complex<double> p_sc = 0., p_sc_add = 0.;

    auto J = [&] () -> std::complex<double> {
                    return ((double)l* sph_bessel(l,k*a) / (k*a) ) - sph_bessel(l+1,k*a);
                };

    auto H = [&] () -> std::complex<double> {
                    return ((double)l* sph_hankel_1(l,k*a) / (k*a) ) - sph_hankel_1(l+1,k*a);
                };

    p_sc_add = - p0 * J()* sph_hankel_1(l,k*r) * legendre_p(l,x[0]/r)/ H();
    p_sc += p_sc_add;

    while ((l<N) && (abs(p_sc_add/p_sc) > tol)) {
        ++l;
        std::complex<double>  A_l = (2*l+1.0)*std::pow(one_i,l)*J()/H();
        p_sc_add = -p0 * A_l * sph_hankel_1(l,k*r) * legendre_p(l,x[0]/r);
        p_sc += p_sc_add;
    }

        values[0] = sin(k*x[0]) + p_sc.imag();
}

};

PYBIND11_MODULE(SIGNATURE, m) {
py::class_<ScatteringExact_Im, std::shared_ptr<ScatteringExact_Im>, dolfin::Expression>
(m, "ScatteringExact_Im")
.def(py::init<>())
.def_readwrite("p0", &ScatteringExact_Im::p0)
.def_readwrite("k", &ScatteringExact_Im::k)
.def_readwrite("a", &ScatteringExact_Im::a)
;
}
