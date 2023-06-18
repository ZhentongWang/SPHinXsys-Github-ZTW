/**
 * @file 	2d_FVM_double_mach_reflection.h
 * @brief 	This is a test to show the double mach reflection in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_DOUBLE_MACH_REFLECTION_H
#define FVM_DOUBLE_MACH_REFLECTION_H
#include "common_shared_FVM_classes.h" // shared eulerian classes for weakly-compressible and compressible fluid in FVM.
#include "common_compressible_FVM_classes.hpp" // eulerian classes for compressible fluid in FVM only.
#include "common_shared_eulerian_classes.h" // shared eulerian classes for weakly-compressible and compressible fluid.
#include "common_compressible_eulerian_classes.h" // eulerian classes for weakly compressible fluid only.
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0;                           /**< Computation domain length. */
Real DH = 1.0;                           /**< Computation domain height. */
Real particle_spacing_ref = 1.0 / 240.0; /**< Initial reference particle spacing. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_one = 1.4;                         /**< initial density of one fluid. */
Real u_one = 0.0;                            /**< initial velocity of one fluid in X axis. */
Real v_one = 0.0;                            /**< initial velocity of one fluid in Y axis. */
Real p_one = 1.0;                            /**< initial pressure of one fluid. */
Real rho0_another = 8.0;                     /**< initial density of another. */
Real u_another = 8.25 * sin(3.14159 / 3.0);  /**< initial velocity of another in X axis. */
Real v_another = -8.25 * cos(3.14159 / 3.0); /**< initial velocity of another in Y axis. */
Real p_another = 140.2 / 1.2;                /**< initial pressure of another. */
Real heat_capacity_ratio = 1.4;              /**< heat capacity ratio. */
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string double_mach_reflection_mesh_fullpath = "./input/double_mach_reflection_0.05.msh";
// 
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> CreatComputationDomian()
{
    //geometry
    std::vector<Vecd> computation_domain;
    computation_domain.push_back(Vecd(0.0, 0.0));
    computation_domain.push_back(Vecd(0.0, DH));
    computation_domain.push_back(Vecd(DL, DH));
    computation_domain.push_back(Vecd(DL, 0.0));
    computation_domain.push_back(Vecd(0.0, 0.0));
    return computation_domain;
}
class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
	{
		MultiPolygon wave_block(CreatComputationDomian());
		add<MultiPolygonShape>(wave_block, "WaveBlock");
	}
};
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class DMFInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit DMFInitialCondition(SPHBody &sph_body)
        : FluidInitialCondition(sph_body), pos_(particles_->pos_), vel_(particles_->vel_),
          rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure"))
    {
        particles_->registerVariable(mom_, "Momentum");
        particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
        particles_->registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
        particles_->registerVariable(E_, "TotalEnergy");
        particles_->registerVariable(dE_dt_, "TotalEnergyChangeRate");
        particles_->registerVariable(dE_dt_prior_, "OtherEnergyChangeRate");
        gamma_ = heat_capacity_ratio;
    };
    void update(size_t index_i, Real dt)
    {
        if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - 1.0 / 6.0))
        {
            /** initial left wave pressure,momentum and energy profile */
            rho_[index_i] = rho0_another;
            p_[index_i] = p_another;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i][0] = u_another;
            vel_[index_i][1] = v_another;
            mom_[index_i] = rho_[index_i] * vel_[index_i];
            E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
        }
        else
        {
            rho_[index_i] = rho0_one;
            p_[index_i] = p_one;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i][0] = u_one;
            vel_[index_i][1] = v_one;
            mom_[index_i] = rho_[index_i] * vel_[index_i];
            E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
        }
    }

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
    StdLargeVec<Real> E_, dE_dt_, dE_dt_prior_;
    Real gamma_;
};

//----------------------------------------------------------------------
//	DMFBoundaryConditionSetup
//----------------------------------------------------------------------
class DMFBoundaryConditionSetup : public LocalDynamics, public DataDelegateInnerInFVM<BaseParticles>
{
public:
	DMFBoundaryConditionSetup(BaseInnerRelationInFVM& inner_relation): 
        LocalDynamics(inner_relation.getSPHBody()), DataDelegateInnerInFVM<BaseParticles>(inner_relation), 
		compressible_fluid_(CompressibleFluid(1.0, 1.4)), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), 
        Vol_(particles_->Vol_), E_(*particles_->getVariableByName<Real>("TotalEnergy")), vel_(particles_->vel_), 
        mom_(*particles_->getVariableByName<Vecd>("Momentum")), pos_(particles_->pos_){};
	virtual ~DMFBoundaryConditionSetup() {};
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
				//rigid wall boundary 
				vel_[index_j] = (vel_[index_i] - e_ij.dot(vel_[index_i])*(e_ij)) + (-e_ij.dot(vel_[index_i])*(e_ij));
				p_[index_j] = p_[index_i];
				rho_[index_j] = rho_[index_i];
				E_[index_j]=E_[index_i];
			}

            if (inner_neighborhood.boundary_type_[n] == 10)
            {
				//given value inlet flow
				Vecd vel_another= Vecd::Zero();
				vel_another[0] = u_another;
				vel_another[1] = v_another;
				Real p_another = 140.2 / 1.2;					/**< initial pressure of another. */
				Real rho_e_another = p_another / (1.4 - 1.0);
				Real E_inlet_another = rho_e_another + 0.5 * rho0_another * vel_another.squaredNorm();

				rho_[index_j] = rho0_another;
				p_[index_j] = p_another;
				vel_[index_j][0] = u_another;
				vel_[index_j][1] = v_another;
				E_[index_j] = E_inlet_another;
            }
                
            if (inner_neighborhood.boundary_type_[n] == 36)
            {
                //Outlet boundary condition
                vel_[index_j] = vel_[index_i];
                p_[index_j] = p_[index_i];
				rho_[index_j] = rho_[index_i];
				E_[index_j]=E_[index_i];

            }

			//Top boundary condition
            if (inner_neighborhood.boundary_type_[n] == 4)
            {
                Real run_time = GlobalStaticVariables::physical_time_;
                Real x_1 = 1.0 / 6.0 + run_time * 10.0 / sin(3.14159 / 3.0);
                if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - x_1))
                {
                    rho_[index_j] = rho0_another;
                    p_[index_j] = p_another;
                    Real rho_e = p_[index_j] / (heat_capacity_ratio - 1.0);
                    vel_[index_j][0] = u_another;
                    vel_[index_j][1] = v_another;
                    E_[index_j] = rho_e + 0.5 * rho_[index_j] * vel_[index_j].squaredNorm();
                }
                else
                {
                    rho_[index_j] = rho0_one;
                    p_[index_j] = p_one;
                    Real rho_e = p_[index_j] / (heat_capacity_ratio - 1.0);
                    vel_[index_j][0] = u_one;
                    vel_[index_j][1] = v_one;
                    E_[index_j] = rho_e + 0.5 * rho_[index_j] * vel_[index_j].squaredNorm();
                }
            }
		}
    };
protected:
	CompressibleFluid compressible_fluid_;
	StdLargeVec<Real>& rho_, & p_, & Vol_, & E_;
	StdLargeVec<Vecd>& vel_, & mom_, & pos_;

};

#endif // FVM_DOUBLE_MACH_REFLECTION_H
