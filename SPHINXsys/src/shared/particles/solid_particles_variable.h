/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	solid_particle_variable.h
* @brief 	Here, we define the algorithm classes for computing derived solid dynamics varabales. 
* @details 	These variable can be added into variable list for state output.   
* @author	Xiangyu Hu
*/

#ifndef SOLID_PARTICLES_VARIABLE_H
#define SOLID_PARTICLES_VARIABLE_H

#include "particle_dynamics_algorithms.h"
#include "body_relation.h"
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
	//----------------------------------------------------------------------
	//		for general solid dynamics variables
	//----------------------------------------------------------------------
	typedef DataDelegateSimple<RealBody, SolidParticles, Solid> SolidDataSimple;

	/**
	 * @class Displacement
	 * @brief computing displacement from current and initial particle position
	 */
	class Displacement : public BaseDerivedVariable<Vecd>, public SolidDataSimple
	{
	public:
		explicit Displacement(SPHBody &sph_body);
		virtual ~Displacement(){};
		void update(size_t index_i, Real dt = 0.0);

	protected:
		StdLargeVec<Vecd> &pos_n_, &pos_0_;
	};

	/**
	 * @class OffsetInitialPosition
	 * @brief offset initial particle position
	 */
	class OffsetInitialPosition : public SolidDataSimple
	{
	public:
		explicit OffsetInitialPosition(SPHBody &sph_body, Vecd &offset);
		virtual ~OffsetInitialPosition(){};
		void update(size_t index_i, Real dt = 0.0);

	protected:
		Vecd offset_;
		StdLargeVec<Vecd> &pos_n_, &pos_0_;
	};

	/**
	 * @class TranslationAndRotation
	 * @brief transfermation on particle position and rotation
	 */
	class TranslationAndRotation : public SolidDataSimple
	{
	public:
		explicit TranslationAndRotation(SPHBody &sph_body, Transformd &transform);
		virtual ~TranslationAndRotation(){};
		void update(size_t index_i, Real dt = 0.0);

	protected:
		Transformd &transform_;
		StdLargeVec<Vecd> &pos_n_, &pos_0_;
	};

	/**
	 * @class NormalDirectionFromBodyShape
	 * @brief normal direction at particles
	 */
	class NormalDirectionFromBodyShape : public SolidDataSimple
	{
	public:
		explicit NormalDirectionFromBodyShape(SPHBody &sph_body);
		virtual ~NormalDirectionFromBodyShape(){};
		void update(size_t index_i, Real dt = 0.0);

	protected:
		Shape &body_shape_;
		StdLargeVec<Vecd> &pos_n_, n_, n_0_;
	};
	
	/**
	 * @class NormalDirectionFromBodyShape
	 * @brief normal direction at particles
	 */
	class NormalDirectionFromShapeAndOp : public SolidDataSimple
	{
	public:
		explicit NormalDirectionFromShapeAndOp(SPHBody &sph_body, const std::string &shape_name);
		virtual ~NormalDirectionFromShapeAndOp(){};
		void update(size_t index_i, Real dt = 0.0);

	protected:
		ShapeAndOp *shape_and_op_;
		Shape *shape_;
		const Real switch_sign_;
		StdLargeVec<Vecd> &pos_n_, n_, n_0_;
	};

	//----------------------------------------------------------------------
	//		for general elastic solid dynamics variables
	//----------------------------------------------------------------------
	typedef DataDelegateSimple<RealBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDataSimple;

	/**
	 * @class VonMisesStress
	 * @brief computing von_Mises_stress
	 */
	class VonMisesStress : public BaseDerivedVariable<Real>, public ElasticSolidDataSimple
	{
	public:
		explicit VonMisesStress(SPHBody &sph_body);
		virtual ~VonMisesStress(){};
		void update(size_t index_i, Real dt = 0.0);

	protected:
		Real rho0_;
		StdLargeVec<Real> &rho_n_;
		StdLargeVec<Matd> &F_, &stress_PK1_;
	};

	/**
	 * @class VonMisesStress
	 * @brief computing von_Mises_stress
	 */
	class VonMisesStrain : public BaseDerivedVariable<Real>, public ElasticSolidDataSimple
	{
	public:
		explicit VonMisesStrain(SPHBody &sph_body);
		virtual ~VonMisesStrain(){};
		void update(size_t index_i, Real dt = 0.0);

	protected:
		StdLargeVec<Matd> &F_;
	};
}
#endif //SOLID_PARTICLES_VARIABLE_H