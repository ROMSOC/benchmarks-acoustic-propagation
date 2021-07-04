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

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <dolfin/function/Expression.h>
#include <iostream>

class exp_uex_re : public dolfin::Expression 
{ public:
    double r1, r2, r3;
    double k_f, rho_f, omega2;
    std::complex<double> k_p, rho_p;
    Eigen::VectorXcd coeffs;

    // Constructor
    exp_uex_re() : dolfin::Expression(3), 
                   r1(0.0), r2(0.0), r3(0.0), omega2(0.0),
                   k_f(0.0), k_p(0.0), rho_f(0.0), rho_p(0.0) {} 

    void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const
    {   auto r = x.norm();
        std::complex<double> one_i(0.0, 1.0);
        std::complex<double> dp_dr, u_mag;

        if (r < r2) { 
        dp_dr = ( coeffs[0] * exp(-one_i * k_f * r) * (-1.0 - one_i * k_f * r) 
                + coeffs[1] * exp( one_i * k_f * r) * (-1.0 + one_i * k_f * r) 
                ) / (r*r);
        u_mag = dp_dr / (rho_f * omega2);
        
        } else if (r < r3) { 
        dp_dr = ( coeffs[2] * exp(-one_i * k_p * r) * (-1.0 - one_i * k_p * r)
                + coeffs[3] * exp( one_i * k_p * r) * (-1.0 + one_i * k_p * r)
                ) / (r*r);
        u_mag = dp_dr / (rho_p * omega2);
        
        } else {
        dp_dr = coeffs[4] * exp(one_i * k_f * r) * (-1.0 + one_i * k_f * r) / (r*r);
        u_mag = dp_dr / (rho_f * omega2);
        }
        
        values[0] = u_mag.real() * x[0] / r; 
        values[1] = u_mag.real() * x[1] / r;
        values[2] = u_mag.real() * x[2] / r;
    }
};

PYBIND11_MODULE(SIGNATURE, m) {
    pybind11::class_<exp_uex_re, std::shared_ptr<exp_uex_re>, dolfin::Expression>
    (m, "exp_uex_re")
    .def(pybind11::init<>())
    .def_readwrite("coeffs", &exp_uex_re::coeffs)
    .def_readwrite("r1", &exp_uex_re::r1)
    .def_readwrite("r2", &exp_uex_re::r2)
    .def_readwrite("r3", &exp_uex_re::r3)
    .def_readwrite("k_f", &exp_uex_re::k_f)
    .def_readwrite("k_p", &exp_uex_re::k_p)
    .def_readwrite("rho_f", &exp_uex_re::rho_f)
    .def_readwrite("rho_p", &exp_uex_re::rho_p)
    .def_readwrite("omega2", &exp_uex_re::omega2)
    ;
}