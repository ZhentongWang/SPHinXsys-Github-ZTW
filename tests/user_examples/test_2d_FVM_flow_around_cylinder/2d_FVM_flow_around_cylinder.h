/**
 * @file 	2d_FVM_flow_around_cylinder.h
 * @brief 	This is a test to show the flow around cylinder case in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_FLOW_AROUND_CYLINDER_H
#define FVM_FLOW_AROUND_CYLINDER_H
#include "common_shared_FVM_classes.h" // shared classes for weakly-compressible and compressible fluid in FVM.
#include "common_weakly_compressible_FVM_classes.hpp" // classes for weakly compressible fluid only in FVM.
#include "common_shared_eulerian_classes.h" // shared eulerian classes for weakly-compressible and compressible fluid.
#include "common_weakly_compressible_eulerian_classes.h" // eulerian classes for weakly compressible fluid only.
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 50.0;						   /**< Channel length. */
Real DH = 30.0;						   /**< Channel height. */
Real resolution_ref = 1.0 / 5.0;	   /**< Initial reference particle spacing. */
Real DL_sponge = 2.0;				   /**< Sponge region to impose inflow condition. */
Real DH_sponge = 2.0;				   /**< Sponge region to impose inflow condition. */
Real cylinder_radius = 1.0;			   /**< Radius of the cylinder. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;										 /**< Density. */
Real U_f = 1.0;											 /**< freestream velocity. */
Real c_f = 10.0 * U_f;									 /**< Speed of sound. */
Real Re = 100.0;										 /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * cylinder_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string zero_three_flow_around_cylinder_mesh_file_fullpath = "./input/fluent_0.3.msh";
// 
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
	water_block_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
	water_block_shape.push_back(Vecd(DL, DH + DH_sponge));
	water_block_shape.push_back(Vecd(DL, -DH_sponge));
	water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

	return water_block_shape;
}
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
	{
		MultiPolygon water_block(createWaterBlockShape());
		add<MultiPolygonShape>(water_block, "WaterBlock");
	}
};

//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class WeaklyCompressibleFluidInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit WeaklyCompressibleFluidInitialCondition(SPHBody &sph_body)
        : FluidInitialCondition(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure"))
    {
        particles_->registerVariable(mom_, "Momentum");
        particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
        particles_->registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
    };
    void update(size_t index_i, Real dt){};

  protected:
    StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
    StdLargeVec<Real> &rho_, &p_;
};

//----------------------------------------------------------------------
//	DMFBoundaryConditionSetup
//----------------------------------------------------------------------
class FACBoundaryConditionSetup : public LocalDynamics, public DataDelegateInnerInFVM<BaseParticles>
{
public:
	FACBoundaryConditionSetup(BaseInnerRelationInFVM& inner_relation): 
        LocalDynamics(inner_relation.getSPHBody()), DataDelegateInnerInFVM<BaseParticles>(inner_relation), 
		compressible_fluid_(CompressibleFluid(1.0, 1.4)), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), 
        Vol_(particles_->Vol_), vel_(particles_->vel_), mom_(*particles_->getVariableByName<Vecd>("Momentum")), pos_(particles_->pos_),
		fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->base_material_)){};
	virtual ~FACBoundaryConditionSetup() {};
	void update(size_t index_i, Real dt = 0.0)
    {
        NeighborhoodInFVM& inner_neighborhood = inner_configuration_in_FVM_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
			Vecd& e_ij = inner_neighborhood.e_ij_[n];

			if (inner_neighborhood.boundary_type_[n] == 3)
			{
				//non-slip wall boundary
				vel_[index_j] = -vel_[index_i];
				p_[index_j] = p_[index_i];
				rho_[index_j] = rho_[index_i];
			}
			if (inner_neighborhood.boundary_type_[n] == 9)
			{
				//Far-field boundary
				Vecd far_field_velocity(1.0, 0.0);
				Real far_field_density = 1.0;
				Real far_field_pressure = fluid_.getPressure(far_field_density);

				vel_[index_j] = far_field_velocity;
				p_[index_j] = far_field_pressure;
				rho_[index_j] = far_field_density;
			}
		}
    };
protected:
	CompressibleFluid compressible_fluid_;
	StdLargeVec<Real>& rho_, & p_, & Vol_;
	StdLargeVec<Vecd>& vel_, & mom_, & pos_;
	Fluid& fluid_;

};

#endif // EULERIAN_FLOW_AROUND_CYLINDER_H
