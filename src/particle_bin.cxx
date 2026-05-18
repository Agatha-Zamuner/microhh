    /*
     		    theta_soil.data(),
    * MicroHH
    * Copyright (c) 2011-2024 Chiel van Heerwaarden
    * Copyright (c) 2011-2024 Thijs Heus
    * Copyright (c) 2014-2024 Bart van Stratum
    *
    * This file is part of MicroHH
    *
    * MicroHH is free software: you can redistribute it and/or modify
    * it under the terms of the GNU General Public License as published by
    * the Free Software Foundation, either version 3 of the License, or
    * (at your option) any later version.

    * MicroHH is distributed in the hope that it will be useful,
    * but WITHOUT ANY WARRANTY; without even the implied warranty of
    * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    * GNU General Public License for more details.

    * You should have received a copy of the GNU General Public License
    * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
    */

    #include <stdexcept>
    #include <iostream>

    #include "master.h"
    #include "input.h"
    #include "grid.h"
    #include "soil_grid.h"
    #include "fields.h"
    #include "constants.h"
    #include "netcdf_interface.h"
    #include "timeloop.h"
    #include "constants.h"
    #include "boundary.h"
    #include "soil_field3d.h"
    #include "particle_bin.h"

    namespace
    {
        template<typename TF>
        void settle_particles(
            TF* const restrict st,
            const TF* const restrict s,
            const TF* const dzhi,
            const TF w_particle,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
        {
            // Simple upwind advection, like in subsidence.
            for (int k=kstart; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jstride + k*kstride;
                st[ijk] -= w_particle * (s[ijk+kstride]-s[ijk])*dzhi[k+1];
            }
        }
        
        template<typename TF>
        void calc_dust_emission(
            TF* const restrict flux_bot,
	    const TF* const restrict theta_soil_top, 
            const TF* const restrict ustar,
            const TF* const restrict TFV_roughness_corr,
            const TF* const restrict eta_f,
            const TF* const restrict eta_c,
            const TF* const restrict sigma_p,
            const TF* const restrict m_ps,
            const TF* const restrict beta_sal,
            const TF* const restrict sand_size,
            const int* const restrict sand_indices,
            const int dust_bin_index,
            const int n_sand,
            const int n_bins,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
        {
            const TF c_y = 5e-5;
            const TF c_q = 2.61;
            const TF rho_a = 1.225;
            const TF rho_bd = 1000.0;
            const TF g = 9.81;
            const TF rho_p = 2650.0;
            const TF p_pcrust = 30000.0;
            const TF k = 0.7;
            const TF theta_r = 0.02; 
            const TF alpha_w = 33.03; 
            const TF beta_w = 0.71; 


            for (int j = jstart; j < jend; ++j)
                for (int i = istart; i < iend; ++i)
            {
                    const int ij = i + j*jstride;
                    const TF u_star = ustar[ij];
                    const TF U = TF(10) * u_star;
                   
	    	    TF f_w; 
		    if (theta_soil_top[ij] > theta_r)
			f_w = std::sqrt(TF(1) + alpha_w * std::pow(theta_soil_top[ij] - theta_r, beta_w));
		    else
			f_w = TF(1);

		    TF TFV_min = TFV_roughness_corr[0] * f_w;
		    for (int b = 1; b < n_bins; ++b)
			TFV_min = std::min(TFV_min, TFV_roughness_corr[b] * f_w);

		    const TF gamma_aggregates = std::exp(-k * std::pow(u_star - TFV_min, TF(3)));

                    TF F_dust_bin = TF(0);

                    for (int s = 0; s < n_sand; ++s)
                    {
                            const int sb = sand_indices[s];
			    const TF TFV_final = TFV_roughness_corr[sb] * f_w; 

                            TF Q = TF(0);
                            if (u_star > TFV_final)
                                Q = c_q * (rho_a/g) * u_star*u_star*u_star * (TF(1) + TFV_final/u_star)* (TF(1) - (TFV_final *TFV_final/(u_star*u_star)));

                            const TF omega = (U*U * sand_size[s]/(beta_sal[s]*beta_sal[s]))* (TF(0.24) + TF(0.21) * U * std::sqrt(rho_p/p_pcrust));

                            const TF term1 = (TF(1) - gamma_aggregates) + (gamma_aggregates * sigma_p[dust_bin_index]);
                            const TF term2 = (g * Q) / (u_star*u_star * m_ps[s]);
                            const TF term3 = rho_bd * eta_f[dust_bin_index] * omega;
                            const TF term4 = eta_c[dust_bin_index] * m_ps[s];

                    //if (i == istart && j == jstart)
                        //std::cout << "d=" << dust_bin_index << "eta_f=" << eta_f[dust_bin_index] << "eta_c=" << eta_c[dust_bin_index] << "\n";
                
                    F_dust_bin += c_y * term1 * term2 * (term3 + term4);

               // if (i == istart && j == jstart)
                      // std::cout << "d=" << dust_bin_index << "F=" << F_dust_bin << "\n";	

                    }
                    flux_bot[ij] += F_dust_bin/rho_a;
            }

        }

    }

    template<typename TF>
    Particle_bin<TF>::Particle_bin(Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soil_gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), soil_grid(soil_gridin), fields(fieldsin)
    {
        sw_particle = inputin.get_item<bool>("particle_bin", "sw_particle", "", false);

        if (sw_particle)
        {
            particle_list = inputin.get_list<std::string>("particle_bin", "particle_list", "", std::vector<std::string>());

            // Read gravitational settling velocities.
            for (auto& scalar : particle_list)
            {
                w_particle.emplace(scalar, inputin.get_item<TF>("particle_bin", "w_particle", scalar));

                // Raise error if any of the velocities is positive.
                if (w_particle.at(scalar) > 0)
                    throw std::runtime_error("Gravitational settling velocities need to be negative!");
            }

            // Constraint on time stepping.
            cfl_max = inputin.get_item<TF>("particle_bin", "cfl_max", "", 1.2);
        }
    }


    template<typename TF>
    Particle_bin<TF>::~Particle_bin()
    {
    }


    template<typename TF>
    void Particle_bin<TF>::init(Netcdf_handle& input_nc)
    {
        if (!sw_particle)
            return;

        // Get dimensions lookup table and resize.
        if (input_nc.group_exists(("particle_bin")))
        {
            Netcdf_handle& nc_group = input_nc.get_group("particle_bin");

            if (nc_group.dimension_exists("dim_x"))
                dim_x = nc_group.get_dimension_size("dim_x");
            if (nc_group.dimension_exists("dim_y"))
                dim_y = nc_group.get_dimension_size("dim_y");
            if (nc_group.dimension_exists("n_bins"))
                n_bins = nc_group.get_dimension_size("n_bins");
            if (nc_group.dimension_exists("n_dust"))
                n_dust = nc_group.get_dimension_size("n_dust");
            if (nc_group.dimension_exists("n_sand"))
                n_sand = nc_group.get_dimension_size("n_sand");

            master.print_message("Particle_bin lookup table shape: (%d, %d)\n", dim_x, dim_y);
        }
    }

    template<typename TF>
    void Particle_bin<TF>::create(Timeloop<TF>& timeloop, Netcdf_handle& input_nc)
    {
        if (!sw_particle)
            return;

        auto& gd = grid.get_grid_data();

        // Calculate fixed maximum time step.
        TF dz_min = TF(Constants::dbig);
        for (int k=gd.kstart; k<gd.kend; ++k)
            dz_min = std::min(dz_min, gd.dz[k]);

        TF w_max = -TF(Constants::dbig);
	for (auto& w : w_particle)
            w_max = std::max(w_max, std::abs(w.second));

        const double dt_max = cfl_max / w_max * dz_min;
        idt_max = convert_to_itime(dt_max);

        // Declare nc_group here so it's available for everything below
        Netcdf_handle& nc_group = input_nc.get_group("particle_bin");

        TFV_roughness_corr.resize(n_bins);
        eta_f.resize(n_dust);
        eta_c.resize(n_dust);
        sigma_p.resize(n_dust);
        m_ps.resize(n_sand);
        beta_sal.resize(n_sand);
        w_terminal.resize(n_bins);
        dust_indices.resize(n_dust);
        sand_indices.resize(n_sand);
        particle_size.resize(n_bins);
        dust_size.resize(n_dust);
        sand_size.resize(n_sand);

        nc_group.get_variable(TFV_roughness_corr, "TFV_roughness_corr", {0}, {n_bins});
        nc_group.get_variable(eta_f, "eta_f", {0}, {n_dust});
        nc_group.get_variable(eta_c, "eta_c", {0}, {n_dust});
        nc_group.get_variable(sigma_p, "sigma_p", {0}, {n_dust});
        nc_group.get_variable(m_ps, "m_ps", {0}, {n_sand});
        nc_group.get_variable(beta_sal, "beta_sal", {0}, {n_sand});
        nc_group.get_variable(w_terminal, "w_terminal", {0}, {n_bins});
        nc_group.get_variable(dust_indices, "dust_indices", {0}, {n_dust});
        nc_group.get_variable(sand_indices, "sand_indices", {0}, {n_sand});
        nc_group.get_variable(particle_size, "particle_size", {0}, {n_bins});
        nc_group.get_variable(dust_size, "dust_size", {0}, {n_dust});
        nc_group.get_variable(sand_size, "sand_size", {0}, {n_sand});


    // for (int s=0; s<n_sand; ++s)
        // master.print_message("sand s=%d, index=%d, TFV_final=%.4f\n",
            // s, sand_indices[s], TFV_final[sand_indices[s]]);

    }

    template<typename TF>
    unsigned long Particle_bin<TF>::get_time_limit()
    {
        if (!sw_particle)
            return Constants::ulhuge;

        return idt_max;
    }


    #ifndef USECUDA
    template<typename TF>
    void Particle_bin<TF>::exec(Boundary<TF>& boundary, Stats<TF>& stats)
    {
        if (!sw_particle)
            return;

        auto& gd = grid.get_grid_data();


         // Gravitational settling of particles.
	for (auto& w : w_particle)
            settle_particles<TF>(
                    fields.st.at(w.first)->fld.data(),
                    fields.sp.at(w.first)->fld.data(),
                    gd.dzhi.data(),
                    w.second,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

        // Surface emissions.
        const std::vector<TF>& ustar = boundary.get_ustar();
	
        TF ustar_min = *std::min_element(ustar.begin(), ustar.end());
        TF ustar_max = *std::max_element(ustar.begin(), ustar.end());
    // master.print_message("ustar min=%.4f, max=%.4f\n", ustar_min, ustar_max); 

	auto& sgd = soil_grid.get_grid_data();
	const TF* theta_soil_top = fields.sps.at("theta")->fld.data()  + (sgd.kend - 1) * (sgd.ncells/sgd.kcells);
	
        for (int d=0; d<dust_indices.size(); ++d)
        {
            const std::string& scalar = particle_list[dust_indices[d]];
            auto& flux_bot = fields.sp.at(scalar)->flux_bot;
            std::fill(flux_bot.begin(), flux_bot.end(), TF(0));
            
            calc_dust_emission(
                    fields.sp.at(scalar)->flux_bot.data(),
                    theta_soil_top,
		    ustar.data(),
                    TFV_roughness_corr.data(),
                    eta_f.data(),
                    eta_c.data(),
                    sigma_p.data(),
                    m_ps.data(),
                    beta_sal.data(),
                    sand_size.data(),
                    sand_indices.data(), 
                    d,
                    n_sand,
                    n_bins,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
    
     	
       
        }
    }


    #endif

    #ifdef FLOAT_SINGLE
    template class Particle_bin<float>;
    #else
    template class Particle_bin<double>;
    #endif


