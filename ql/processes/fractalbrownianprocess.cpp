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

#include <ql/processes/fractalbrownianprocess.hpp>
#include <ql/processes/eulerdiscretization.hpp>

namespace QuantLib {

	FractalBrownianMotionProcess::FractalBrownianMotionProcess(
														double initialValue,
                                                        double mue,
                                                        double sigma,
														double hurst,
														int length)
    : StochasticProcess1D(boost::shared_ptr<discretization>(
                                                    new EulerDiscretization)),
      initialValue_(initialValue), mue_(mue), sigma_(sigma), hurst_(hurst), length_(length) 
	{
	// Plan, when a fractal process is instanziated, a fgn and a fbm with length length_ and Hurst-Exponent hurst_ 
	// should be generated immediately. Later evolve will just return the values of the process previously computed.
	// In a first attempt, a fixed RNG should be used.

		// Define RN generator (Gaussian via Box-Muller) and sequence generator:
		BoxMullerGaussianRng<MersenneTwisterUniformRng> gaussianRng(MersenneTwisterUniformRng(0));
		RandomSequenceGenerator<BoxMullerGaussianRng<MersenneTwisterUniformRng>> gaussianSequence(length_, gaussianRng);

		// Get uncorrelated sequence:
		seq_type gn = gaussianSequence.nextSequence();

		// TODO: Implement FGN-Class and correlated sequence...

		seq_type fgn = gn, fbm = gn;

		// Cumulatively compute fBM from fGn:
		for (int i = 1; i < length; i++)
		{
			fbm.value[i] = fbm.value[i - 1] + fgn.value[i];
		}
	}

    Real FractalBrownianMotionProcess::x0() const {
        return initialValue_;
    }

    Real FractalBrownianMotionProcess::drift(Time, Real x) const {
        return mue_ * x;
    }

    Real FractalBrownianMotionProcess::diffusion(Time, Real x) const {
        return sigma_ * x;
    }

	Real FractalBrownianMotionProcess::hurst() const {
		return hurst_;
	}

	Size FractalBrownianMotionProcess::length() const {
		return length_;
	}

}
