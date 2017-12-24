/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017, Gregor Wergen

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file fractalbrownianprocess.hpp
    \brief Fractal Brownian-motion process
*/

#ifndef quantlib_fractal_brownian_process_hpp
#define quantlib_fractal_brownian_process_hpp

#include <ql/stochasticprocess.hpp>
#include <ql/math/randomnumbers/boxmullergaussianrng.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include <vector>
#include <ql/math/randomnumbers/randomsequencegenerator.hpp>

namespace QuantLib {

    //! Fractal brownian-motion process

    class FractalBrownianMotion {
      public:
        FractalBrownianMotion(double hurst);

		Real x0() const;
		Real drift(Time t, Real x) const;
		Real diffusion(Time t, Real x) const;
		Real hurst() const;
		Size length() const;
		typedef Sample<Real> sample_type;
		typedef Sample<std::vector<Real> > seq_type;
      protected:
		  double hurst_;
    };

}


#endif
