/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*                                                                         *
 * ------------------------------------------------------------------------*/
#pragma once

#include "common_compressible_FVM_classes.h"


namespace SPH
{
	//=================================================================================================//
	template <class RiemannSolverType>
	BaseIntegration1stHalfInFVM<RiemannSolverType>::BaseIntegration1stHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter)
    : BaseIntegrationInCompressibleFVM(inner_relation), riemann_solver_(compressible_fluid_, compressible_fluid_, limiter_parameter) {}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration1stHalfInFVM<RiemannSolverType>::interaction(size_t index_i, Real dt)
	{
		CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], E_[index_i]);
		Vecd momentum_change_rate = dmom_dt_prior_[index_i];
		NeighborhoodInFVM& inner_neighborhood = inner_configuration_in_FVM_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
			Vecd& e_ij = inner_neighborhood.e_ij_[n];

			CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], E_[index_j]);
			CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

			momentum_change_rate -= 2.0 * dW_ijV_j *
			((interface_state.rho_ * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij;
		}
		dmom_dt_[index_i] = momentum_change_rate;
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration1stHalfInFVM<RiemannSolverType>::update(size_t index_i, Real dt)
	{
		mom_[index_i] += dmom_dt_[index_i] * dt;
		vel_[index_i] = mom_[index_i] / rho_[index_i];
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	BaseIntegration2ndHalfInFVM<RiemannSolverType>::BaseIntegration2ndHalfInFVM(BaseInnerRelationInFVM &inner_relation, Real limiter_parameter)
    : BaseIntegrationInCompressibleFVM(inner_relation), riemann_solver_(compressible_fluid_, compressible_fluid_,limiter_parameter) {}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration2ndHalfInFVM<RiemannSolverType>::interaction(size_t index_i, Real dt)
	{
		CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], E_[index_i]);
		Real density_change_rate = 0.0;
		Real energy_change_rate = dE_dt_prior_[index_i];
		NeighborhoodInFVM& inner_neighborhood = inner_configuration_in_FVM_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Vecd& e_ij = inner_neighborhood.e_ij_[n];
			Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

			CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], E_[index_j]);
			CompressibleFluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

			density_change_rate -= 2.0 * dW_ijV_j * (interface_state.rho_ * interface_state.vel_).dot(e_ij);
			energy_change_rate -= 2.0 * dW_ijV_j * (interface_state.E_ * interface_state.vel_ + interface_state.p_ * interface_state.vel_).dot(e_ij);
		}
		drho_dt_[index_i] = density_change_rate;
		dE_dt_[index_i] = energy_change_rate;
	}
	//=================================================================================================//
	template <class RiemannSolverType>
	void BaseIntegration2ndHalfInFVM<RiemannSolverType>::update(size_t index_i, Real dt)
	{
		E_[index_i] += dE_dt_[index_i] * dt ;
		rho_[index_i] += drho_dt_[index_i] * dt;
		Real rho_e = E_[index_i] - 0.5 * mom_[index_i].squaredNorm() / rho_[index_i];
		p_[index_i] = compressible_fluid_.getPressure(rho_[index_i], rho_e);
	}
	//=================================================================================================//
}
//=================================================================================================//