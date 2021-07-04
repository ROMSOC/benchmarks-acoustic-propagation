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

class exp_pex_im : public dolfin::Expression 
{ public:
    double r1, r2, r3;
    double k_f;
    std::complex<double> k_p;
    Eigen::VectorXcd coeffs;

    exp_pex_im() : dolfin::Expression(), r1(0.0), r2(0.0), r3(0.0), k_f(0.0), k_p(0.0) {}

    void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const 
    {  auto r = x.norm();
        std::complex<double> one_i(0.0, 1.0);
        std::complex<double> p; 
        
        if (r < r2) p = (coeffs[0] * exp(-one_i * k_f * r) + coeffs[1] * exp(one_i * k_f * r)) / r;
        else if (r < r3) p = (coeffs[2] * exp(-one_i * k_p * r) + coeffs[3] * exp(one_i * k_p * r)) / r;
        else p = coeffs[4] * exp(one_i * k_f * r) / r;

        values[0] = p.imag();
    }
};

PYBIND11_MODULE(SIGNATURE, m) {
    pybind11::class_<exp_pex_im, std::shared_ptr<exp_pex_im>, dolfin::Expression>
    (m, "exp_pex_im")
    .def(pybind11::init<>())
    .def_readwrite("coeffs", &exp_pex_im::coeffs)
    .def_readwrite("r1", &exp_pex_im::r1)
    .def_readwrite("r2", &exp_pex_im::r2)
    .def_readwrite("r3", &exp_pex_im::r3)
    .def_readwrite("k_f", &exp_pex_im::k_f)
    .def_readwrite("k_p", &exp_pex_im::k_p)
    ;
}
