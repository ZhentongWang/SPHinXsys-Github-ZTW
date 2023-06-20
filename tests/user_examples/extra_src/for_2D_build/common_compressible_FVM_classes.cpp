
#include "common_compressible_FVM_classes.h"

namespace SPH
{
	//=================================================================================================//
	CompressibleAcousticTimeStepSizeInFVM::CompressibleAcousticTimeStepSizeInFVM(SPHBody& sph_body)
		: AcousticTimeStepSize(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), vel_(particles_->vel_), 
		compressible_fluid_(CompressibleFluid(1.0, 1.4)) {};
	//=================================================================================================//
	Real CompressibleAcousticTimeStepSizeInFVM::reduce(size_t index_i, Real dt)
	{
		return compressible_fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
	}
	//=================================================================================================//
	Real CompressibleAcousticTimeStepSizeInFVM::outputResult(Real reduced_value)
	{
		return 0.0005 / (reduced_value + TinyReal);
	}
	//=================================================================================================//
}
//=================================================================================================//