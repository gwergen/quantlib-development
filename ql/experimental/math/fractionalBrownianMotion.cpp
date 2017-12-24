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

#include <ql/experimental/math/fractionalBrownianMotion.hpp>

namespace QuantLib {

	FractalBrownianMotion::FractalBrownianMotion(
			Real hurst,
			Real length, 
			int n,
			std::string method="cholesky") :
		hurst_(hurst), length_(length), increments_(std::vector<Real>(n), 0.0), process_(std::vector<Real>(n), 0.0), method_(method), n_(n)
	{
		QL_REQUIRE(method_ == "cholesky" || method_ == "hosking" || method_ == "daviesharte",
			"method " << method_ << " not supported, must be 'hosking', 'cholesky' or 'daviesharte'!");
		QL_REQUIRE(hurst_ > 0.0 && hurst_ < 1.0,
			"FBM for a Hurst-Exponent of  " << hurst_ << " is not defined, it must be > 0.0 and < 1.0!");

		// Define RN generator (Gaussian via Box-Muller) and sequence generator:
		BoxMullerGaussianRng<MersenneTwisterUniformRng> gaussianRng(MersenneTwisterUniformRng(0));
		RandomSequenceGenerator<BoxMullerGaussianRng<MersenneTwisterUniformRng>> gaussianSequence(n_, gaussianRng);

		// Get uncorrelated sequence:
		seq_type gn = gaussianSequence.nextSequence();
		seq_type fgn = gn, fbm = gn;

		if (method_ == "cholesky") fgn = cholesky(gn);
		else if (method_ == "hosking") fgn = hosking(gn);
		else fgn = davieshartie(gn);

		// Cumulatively compute fBM from fGn:
		for (int i = 1; i < n_; i++) fbm.value[i] = fbm.value[i - 1] + fgn.value[i];

		increments_ = fgn;
		process_ = fbm;
	}

	Real FractalBrownianMotion::hurst() const {
		return hurst_;
	}

	Real FractalBrownianMotion::length() const {
		return length_;
	}

	int FractalBrownianMotion::steps() const {
		return n_;
	}

	FractalBrownianMotion::seq_type FractalBrownianMotion::getIncrements() {
		return increments_;
	}

	FractalBrownianMotion::seq_type FractalBrownianMotion::getProcess() {
		return process_;
	}

	FractalBrownianMotion::seq_type FractalBrownianMotion::cholesky(seq_type gn) {
		Real increment = length_ / n_;
		Real scale = pow(increment, hurst_);
		if(hurst_ != 0.5) { 
			seq_type fgn=gn;
			for (int i = 0; i < n_;i++) fgn.value[i] = 0.0;
			Matrix G(n_, n_);
			for (int i = 0; i < n_;i++) {
				for (int j = 0; j < n_;j++) {
					G(i, j) = autocovariance(i - j);
				}
			}
			Matrix C = CholeskyDecomposition(G);
			for (int i = 0; i < n_;i++) fgn.value[i] = 0;
			for (int i = 0; i < n_;i++) for (int j = 0; j < n_;j++) {
				fgn.value[i] += C(i,j) * gn.value[i];
			}
			for (int i= 0; i<n_; i++) fgn.value[i] *= scale;
			return fgn;
		}
		return gn;
	}

	FractalBrownianMotion::seq_type FractalBrownianMotion::hosking(seq_type gn) {
		Real increment = length_ / n_;
		Real scale = pow(increment, hurst_);
		if (hurst_ != 0.5) {
			seq_type fgn(std::vector<Real>(n_, 0.0), 1.0);
			std::vector<Real> phi(n_, 0.0);
			std::vector<Real> psi(n_, 0.0);
			std::vector<Real> cov(n_, 0.0);
			for (int i=0; i < n_; i++) cov[i] = autocovariance(i);
			double v = 1.0;
			fgn.value[0] = gn.value[0];

			for (int i=1; i < n_; i++) {
				phi[i - 1] = cov[i];
				for (int j = 0; j < i - 1; j++) {
					psi[j] = phi[j];
					phi[i - 1] -= psi[j] * cov[i - j - 1];
				}
				phi[i - 1] /= v;
				for (int j = 0; j < i - 1; j++) {
					phi[j] = psi[j] - phi[i - 1] * psi[i - j - 2];
				}
				v *= (1 - phi[i - 1] * phi[i - 1]);
				for (int j = 0; j < i; j++) {
					fgn.value[i] += phi[j] * fgn.value[i - j - 1];
				}
				fgn.value[i] += sqrt(v) * gn.value[i];
			}
			for (int i = 0; i<n_; i++) fgn.value[i] *= scale;
			return fgn;
		}
		for (int i = 0; i<n_; i++) gn.value[i] *= scale;
		return gn;
	}

	FractalBrownianMotion::seq_type FractalBrownianMotion::davieshartie(seq_type gn) {
		Real increment = length_ / n_;
		Real scale = pow(increment, hurst_);
		if (hurst_ != 0.5) {
			seq_type fgn = gn;

			std::vector<Real> row_component(n_ - 1, 0.0);
			std::vector<Real> reverse_component(n_ - 1, 0.0);
			std::vector<Real> row(2*n_, 0.0);

			for (int i = 1;i < n_;i++) row_component[i - 1] = autocovariance(i);
			for (int i = 1;i < n_;i++) reverse_component[i - 1] = autocovariance(-i);

			row[0] = autocovariance(0);
			for (int i = 1;i < n_;i++) row[i] = row_component[i-1];
			row[n_] = 0.0;
			for (int i = n_ + 1; i < 2 * n_;i++) row[i] = reverse_component[i - n_ - 1];

			/*
			
			Remains to be implemented...
			
			*/

			return fgn;
		}
		return gn;
	}

	Real FractalBrownianMotion::autocovariance(int distance) {
		return 0.5 * (pow(abs(distance - 1), 2 * hurst_)
			- 2 * pow(abs(distance), 2 * hurst_)
			+ pow(abs(distance + 1), 2 * hurst_));
	}
}