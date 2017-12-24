/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 Ferdinando Ametrano
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006 StatPro Italia srl

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

/*! \file fractalpathgenerator.hpp
    \brief Generates random paths using a sequence generator
*/

#ifndef quantlib_montecarlo_fractal_path_generator_hpp
#define quantlib_montecarlo_fractal_path_generator_hpp

#include <ql/methods/montecarlo/brownianbridge.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/experimental/math/fractionalBrownianMotion.hpp>

namespace QuantLib {
    class StochasticProcess;
    class StochasticProcess1D;
    //! Generates random paths using a sequence generator
    /*! Generates random paths with drift(S,t) and variance(S,t)
        using a gaussian sequence generator

        \ingroup mcarlo

        \test the generated paths are checked against cached results
    */
    template <class GSG>
    class FractalPathGenerator {
      public:
        typedef Sample<Path> sample_type;
        // constructors
		FractalPathGenerator(const boost::shared_ptr<StochasticProcess>&,
                      Time length,
                      Size timeSteps,
					  Real hurst,
					  std::string method,
                      const GSG& generator);
		FractalPathGenerator(const boost::shared_ptr<StochasticProcess>&,
                      const TimeGrid& timeGrid,
					  Real hurst,
					  std::string method,
                      const GSG& generator);
        //! \name inspectors
        //@{
        const sample_type& next() const;
        Size size() const { return dimension_; }
        const TimeGrid& timeGrid() const { return timeGrid_; }
		Real hurst() const { return hurst_; }
		std::string method() const { return method_; }
        //@}
      private:
        //const sample_type& next(bool antithetic) const;
        GSG generator_;
        Size dimension_;
        TimeGrid timeGrid_;
        boost::shared_ptr<StochasticProcess1D> process_;
        mutable sample_type next_;
        mutable std::vector<Real> temp_;
		Real hurst_;
		std::string method_;
    };

    // template definitions

    template <class GSG>
	FractalPathGenerator<GSG>::FractalPathGenerator(
                          const boost::shared_ptr<StochasticProcess>& process,
                          Time length,
                          Size timeSteps,
						  Real hurst,
						  std::string method,
                          const GSG& generator)
    : generator_(generator), dimension_(generator_.dimension()), timeGrid_(length, timeSteps), 
	  hurst_(hurst), method_(method), process_(boost::dynamic_pointer_cast<StochasticProcess1D>(process)),
      next_(Path(timeGrid_),1.0), temp_(dimension_) {
        QL_REQUIRE(dimension_==timeSteps,
                   "sequence generator dimensionality (" << dimension_
                   << ") != timeSteps (" << timeSteps << ")");
    }

    template <class GSG>
	FractalPathGenerator<GSG>::FractalPathGenerator(
                          const boost::shared_ptr<StochasticProcess>& process,
                          const TimeGrid& timeGrid,
						  Real hurst,
						  std::string method,
                          const GSG& generator)
    : generator_(generator),
      dimension_(generator_.dimension()), timeGrid_(timeGrid), hurst_(hurst), method_(method),
      process_(boost::dynamic_pointer_cast<StochasticProcess1D>(process)),
      next_(Path(timeGrid_),1.0), temp_(dimension_) {
        QL_REQUIRE(dimension_==timeGrid_.size()-1,
                   "sequence generator dimensionality (" << dimension_
                   << ") != timeSteps (" << timeGrid_.size()-1 << ")");
    }

    template <class GSG>
    const typename FractalPathGenerator<GSG>::sample_type&
    FractalPathGenerator<GSG>::next() const {

        typedef typename GSG::sample_type sequence_type;

		FractalBrownianMotion FBM(hurst_, timeGrid_.back(), timeGrid_.size(), method_);
		Sample<std::vector<Real> > fbn = FBM.getIncrements();
		Path& path = next_.value;

		// Currently, the new fractal path adopts x0 and standard deviation from process_
		path.front() = process_->x0();
		for (Size i = 1; i<path.length(); i++) {
			Time t = timeGrid_[i - 1];
			Time dt = timeGrid_.dt(i - 1);
			Real sigma = process_->stdDeviation(t, path[i - 1], dt);
			Real mu = process_->drift(t, path[i - 1]);
			path[i] = path[i-1] + (fbn.value[i-1]);
		}

		return next_;
    }

}


#endif
