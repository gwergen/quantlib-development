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

/*! \file fractionalBrownianMotion.hpp
\brief Generates fractional brownian motion
*/

#ifndef quantlib_fractional_brownian_motion_hpp
#define quantlib_fractional_brownian_motion_hpp

#include <ql/stochasticprocess.hpp>
#include <ql/math/randomnumbers/boxmullergaussianrng.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include <vector>
#include <math.h>
#include <ql/math/randomnumbers/randomsequencegenerator.hpp>
#include <ql/math/matrix.hpp>
#include <ql/math/matrixutilities/choleskydecomposition.hpp>

namespace QuantLib {

	//! Fractal brownian-motion

	class FractalBrownianMotion {
	public:
		FractalBrownianMotion(Real hurst, Real length, int n, std::string method);

		Real hurst() const;
		Real length() const;
		int steps() const;

		typedef Sample<Real> sample_type;
		typedef Sample<std::vector<Real> > seq_type;

		seq_type getIncrements();
		seq_type getProcess();

	private:
		seq_type cholesky(seq_type seq);
		seq_type hosking(seq_type seq);
		seq_type davieshartie(seq_type seq);
		Real autocovariance(int i);

	protected:
		Real hurst_;
		Real length_;
		int n_;
		seq_type increments_;
		seq_type process_;
		std::string method_;
	};

}

#endif