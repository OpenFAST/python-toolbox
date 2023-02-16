import numpy as np
import os

class AMRWindSimulation:
    '''
    This class is used to help prepare sampling planes for an AMR-Wind
      simulation. The sampling planes will be used to generate inflow
      data for FAST.Farm simulations.
    Specifically, this class contains info from the AMR-Wind input file,
      and it carries out simple calculations about the AMR-Wind
      simulation
    For reference, see https://openfast.readthedocs.io/en/dev/source/user/fast.farm/ModelGuidance.html
    '''

    def __init__(self, wts:dict,
                    dt: float, prob_lo: tuple, prob_hi: tuple, 
                    n_cell: tuple, max_level: int, 
                    incflo_velocity_hh: tuple, 
                    postproc_name='sampling',
                    buffer_lr = [3,6,3,3,2],
                    buffer_hr = 0.6,
                    ds_hr = None, ds_lr = None):
        '''
        Values from the AMR-Wind input file
        Inputs:
          * dt: this should be a fixed dt value from the LES run
          * incflo_velocity_hh: velocity vector, specifically at hub height
          * buffer_lr: buffer for [xmin, xmax, ymin, ymax, zmax] in low-res box, in D
          * buffer_hr: buffer for all directions (constant) in high-res box, in D
        '''
        # Process inputs
        self.wts                = wts
        self.dt                 = dt
        self.prob_lo            = prob_lo
        self.prob_hi            = prob_hi
        self.n_cell             = n_cell
        self.max_level          = max_level
        self.incflo_velocity_hh = incflo_velocity_hh
        self.postproc_name_lr   = f"{postproc_name}_lr"
        self.postproc_name_hr   = f"{postproc_name}_hr"
        self.buffer_lr          = buffer_lr
        self.buffer_hr          = buffer_hr
        self.ds_hr              = ds_hr
        self.ds_lr              = ds_lr

        # Placeholder variables, to be calculated by FFCaseCreation
        self.output_frequency_lr = None
        self.output_frequency_hr = None
        self.sampling_labels_lr = None
        self.sampling_labels_hr = None
        self.nx_lr = None
        self.ny_lr = None
        self.nz_lr = None
        self.xlow_lr = None
        self.xhigh_lr = None
        self.ylow_lr = None
        self.yhigh_lr = None
        self.zlow_lr = None
        self.zhigh_lr = None
        self.zoffsets_lr = None
        self.hr_domains = None

        # Run extra functions
        self._checkInputs()
        self._calc_simple_params()
        self._calc_sampling_params()

    def _checkInputs(self):
        '''
        Check that the AMR-Wind inputs make sense
        '''
        if len(self.prob_lo) != 3:
            raise ValueError(f"prob_lo must contain 3 elements, but it has {len(self.prob_lo)}")
        if len(self.prob_hi) != 3:
            raise ValueError(f"prob_hi must contain 3 elements, but it has {len(self.prob_hi)}")
        if len(self.incflo_velocity_hh) != 3:
            raise ValueError(f"incflo_velocity_hh must contain 3 elements, but it has {len(self.incflo_velocity_hh)}")
        if (self.prob_lo[0] >= self.prob_hi[0]):
            raise ValueError("x-component of prob_lo larger than x-component of prob_hi")
        if (self.prob_lo[1] >= self.prob_hi[1]):
            raise ValueError("y-component of prob_lo larger than y-component of prob_hi")
        if (self.prob_lo[2] >= self.prob_hi[2]):
            raise ValueError("z-component of prob_lo larger than z-component of prob_hi")

    def _calc_simple_params(self):
        '''
        Calculate simulation parameters, given only AMR-Wind inputs
        '''
        # Grid resolution at Level 0
        self.dx0 = (self.prob_hi[0] - self.prob_lo[0]) / self.n_cell[0]
        self.dy0 = (self.prob_hi[1] - self.prob_lo[1]) / self.n_cell[1]
        self.dz0 = (self.prob_hi[2] - self.prob_lo[2]) / self.n_cell[2]
        self.ds0_max = max(self.dx0, self.dy0, self.dz0)

        # Grid resolution at finest refinement level
        self.dx_refine = self.dx0/(2**self.max_level)
        self.dy_refine = self.dy0/(2**self.max_level)
        self.dz_refine = self.dz0/(2**self.max_level)
        self.ds_refine_max = max(self.dx_refine, self.dy_refine, self.dz_refine)

        # Hub height wind speed
        self.vhub = np.sqrt(self.incflo_velocity_hh[0]**2 + self.incflo_velocity_hh[1]**2)

    def _calc_sampling_params(self):
        '''
        Calculate parameters for sampling planes
        '''

        self._calc_sampling_labels()
        self._calc_sampling_time()
        self._calc_grid_resolution()
        self._calc_grid_placement()

    def _calc_sampling_labels(self):
        '''
        Calculate labels for AMR-Wind sampling
        '''
        sampling_labels_lr = ["Low"]
        self.sampling_labels_lr = sampling_labels_lr

        sampling_labels_hr = []
        for turbkey in self.wts:
            if 'name' in self.wts[turbkey].keys():
                wt_name = self.wts[turbkey]['name']
            else:
                wt_name = f'T{turbkey}'
            sampling_labels_hr.append(f"High{wt_name}_inflow0deg")
        
        self.sampling_labels_hr = sampling_labels_hr

    def _calc_sampling_time(self):
        '''
        Calculate timestep values and AMR-Wind plane sampling frequency
        '''

        ## High resolution domain, dt_high_les
        fmax_max = 0
        for turbkey in self.wts:
            fmax_max = max(0, self.wts[turbkey]['fmax'])
        dt_hr_max = 1 / (2 * fmax_max)
        self.dt_high_les = self.dt * np.floor(dt_hr_max/self.dt)  # Ensure that dt_hr is a multiple of the AMR-Wind timestep

        if self.dt_high_les < self.dt:
            raise ValueError(f"AMR-Wind timestep {self.dt} too coarse for high resolution domain! AMR-Wind timestep must be at least {self.dt_high_les} sec.")

        ## Low resolution domain, dt_low_les
        cmeander_min = float("inf")
        Dwake_min = float("inf")
        for turbkey in self.wts:
            cmeander_min = min(cmeander_min, self.wts[turbkey]['Cmeander'])
            Dwake_min = min(Dwake_min, self.wts[turbkey]['D'])  # Approximate D_wake as D_rotor

        dt_lr_max = cmeander_min * Dwake_min / (10 * self.vhub)
        self.dt_low_les = self.dt_high_les * np.floor(dt_lr_max/self.dt_high_les)  # Ensure that dt_lr is a multiple of the high res sampling timestep

        if self.dt_low_les < self.dt:
            raise ValueError(f"AMR-Wind timestep {self.dt} too coarse for low resolution domain! AMR-Wind timestep must be at least {self.dt_low_les} sec.")
        if self.dt_high_les > self.dt_low_les:
            raise ValueError(f"Low resolution timestep ({self.dt_low_les}) is finer than high resolution timestep ({self.dt_high_les})!")

        ## Sampling frequency
        self.output_frequency_hr = int(np.floor(self.dt_high_les/self.dt))
        self.output_frequency_lr = self.output_frequency_hr * np.floor(self.dt_low_les/self.dt_high_les)

        if self.output_frequency_lr % self.output_frequency_hr != 0:
            raise ValueError(f"Low resolution output frequency of {self.output_frequency_lr} not a multiple of the high resolution frequency {self.output_frequency_hr}!")

    def _calc_grid_resolution(self):
        '''
        Calculate sampling grid resolutions
        '''

        ## High resolution domain, ds_hr
        #    ASSUME: FAST.Farm HR zone lies within the region of maxmum AMR-Wind grid refinement
        #    NOTE: ds_hr is calculated independent of any x/y/z requirements,
        #            just blade chord length requirements
        cmax_min = float("inf")
        for turbkey in self.wts:
            cmax_min = min(cmax_min, self.wts[turbkey]['cmax'])
        ds_hr_max = cmax_min

        if self.ds_hr is None:  # Calculate ds_hr if it is not specified as an input
            if ds_hr_max < self.ds_refine_max:
                raise ValueError(f"AMR-Wind grid spacing of {self.ds_refine_max} is too coarse for high resolution domain! The high-resolution domain requires "\
                                 f"AMR-Wind grid spacing to be at least {ds_hr_max} m. If a coarser high-res domain is acceptable, then manually specify the "\
                                 f"high-resolution grid spacing to be at least {self.ds_refine_max} with ds_hr = {self.ds_refine_max}.")
            self.ds_high_les = self.ds_refine_max * np.floor(ds_hr_max/self.ds_refine_max)  # Ensure that ds_hr is a multiple of the refined AMR-Wind grid spacing
            self.ds_hr = self.ds_high_les
        else:
            self.ds_high_les = self.ds_hr

        if self.ds_high_les < self.ds_refine_max:
            raise ValueError(f"AMR-Wind fine grid spacing {self.ds_refine_max} too coarse for high resolution domain! AMR-Wind grid spacing must be at least {self.ds_high_les} m.")

        ## Low resolution domain, ds_lr (s = x/y/z)
        #    ASSUME: FAST.Farm LR zone uses Level 0 AMR-Wind grid spacing
        #    NOTE: ds_lr is calculated independent of any x/y/z requirements,
        #            just time step and velocity requiements
        if self.ds_lr is None:
            ds_lr_max = self.dt_low_les * self.vhub**2 / 15
            self.ds_low_les = self.ds0_max * np.floor(ds_lr_max/self.ds0_max)  # Ensure that ds_lr is a multiple of the coarse AMR-Wind grid spacing
            self.ds_lr = self.ds_low_les
        else:
            self.ds_low_les = self.ds_lr

        if self.ds_low_les < self.ds0_max:
            raise ValueError(f"AMR-Wind coarse grid spacing {self.ds0_max} too coarse for high resolution domain! AMR-Wind grid spacing must be at least {self.ds_low_les} m.")
        if self.ds_low_les % self.ds_high_les != 0:
            raise ValueError(f"Low resolution grid spacing of {self.ds_low_les} not a multiple of the high resolution grid spacing {self.ds_high_les}!")

    def _calc_grid_placement(self):
        '''
        Calculate placement of sampling grids
        '''

        self._calc_hr_grid_placement()
        self._calc_lr_grid_placement()

    def _calc_hr_grid_placement(self):
        '''
        Calculate placement of high resolution grids
        '''
        ### ~~~~~~~~~ Calculate high resolution grid placement ~~~~~~~~~
        hr_domains = {} 
        for turbkey in self.wts:
            wt_x = self.wts[turbkey]['x']
            wt_y = self.wts[turbkey]['y']
            wt_z = self.wts[turbkey]['zhub']
            wt_D = self.wts[turbkey]['D']

            # Calculate minimum/maximum HR domain extents
            xlow_hr_min  = wt_x - self.buffer_hr * wt_D 
            xhigh_hr_max = wt_x + self.buffer_hr * wt_D 
            ylow_hr_min  = wt_y - self.buffer_hr * wt_D 
            yhigh_hr_max = wt_y + self.buffer_hr * wt_D 
            zhigh_hr_max = wt_z + self.buffer_hr * wt_D 

            # Calculate the minimum/maximum HR domain coordinate lengths & number of grid cells
            xdist_hr_min = xhigh_hr_max - xlow_hr_min  # Minumum possible length of x-extent of HR domain
            xdist_hr = self.ds_high_les * np.ceil(xdist_hr_min/self.ds_high_les)  
            nx_hr = int(xdist_hr/self.ds_high_les) + 1

            ydist_hr_min = yhigh_hr_max - ylow_hr_min
            ydist_hr = self.ds_high_les * np.ceil(ydist_hr_min/self.ds_high_les)
            ny_hr = int(ydist_hr/self.ds_high_les) + 1

            zdist_hr = self.ds_high_les * np.ceil(zhigh_hr_max/self.ds_high_les)
            nz_hr = int(zdist_hr/self.ds_high_les) + 1

            # Calculate actual HR domain extent
            #  NOTE: Sampling planes should measure at AMR-Wind cell centers, not cell edges
            xlow_hr = self.ds_high_les * np.floor(xlow_hr_min/self.ds_high_les) - 0.5*self.dx_refine + self.prob_lo[0]%self.ds_high_les
            xhigh_hr = xlow_hr + xdist_hr
            ylow_hr = self.ds_high_les * np.floor(ylow_hr_min/self.ds_high_les) - 0.5*self.dy_refine + self.prob_lo[1]%self.ds_high_les
            yhigh_hr = ylow_hr + ydist_hr
            zlow_hr = 0.5 * self.dz0 / (2**self.max_level)
            zhigh_hr = zlow_hr + zdist_hr
            zoffsets_hr = np.arange(zlow_hr, zhigh_hr+self.ds_high_les, self.ds_high_les) - zlow_hr

            # Check domain extents
            if xhigh_hr > self.prob_hi[0]:
                raise ValueError(f"HR domain point {xhigh_hr} extends beyond maximum AMR-Wind x-extent!")
            if xlow_hr < self.prob_lo[0]:
                raise ValueError(f"HR domain point {xlow_hr} extends beyond minimum AMR-Wind x-extent!")
            if yhigh_hr > self.prob_hi[1]:
                raise ValueError(f"HR domain point {yhigh_hr} extends beyond maximum AMR-Wind y-extent!")
            if ylow_hr < self.prob_lo[1]:
                raise ValueError(f"HR domain point {ylow_hr} extends beyond minimum AMR-Wind y-extent!")
            if zhigh_hr > self.prob_hi[2]:
                raise ValueError(f"HR domain point {zhigh_hr} extends beyond maximum AMR-Wind z-extent!")
            if zlow_hr < self.prob_lo[2]:
                raise ValueError(f"HR domain point {zlow_hr} extends beyond minimum AMR-Wind z-extent!")

            # Save out info for FFCaseCreation
            self.extent_high = self.buffer_hr

            hr_turb_info = {'nx_hr': nx_hr, 'ny_hr': ny_hr, 'nz_hr': nz_hr,
                            'xdist_hr': xdist_hr, 'ydist_hr': ydist_hr, 'zdist_hr': zdist_hr,
                            'xlow_hr': xlow_hr, 'ylow_hr': ylow_hr, 'zlow_hr': zlow_hr,
                            'xhigh_hr': xhigh_hr, 'yhigh_hr': yhigh_hr, 'zhigh_hr': zhigh_hr,
                            'zoffsets_hr': zoffsets_hr}
            hr_domains[turbkey] = hr_turb_info
        self.hr_domains = hr_domains

    def _calc_lr_grid_placement(self):
        '''
        Calculate placement of low resolution grid
        '''

        ### ~~~~~~~~~ Calculate low resolution grid placement ~~~~~~~~~ 
        # Calculate minimum/maximum LR domain extents
        wt_all_x_min = float("inf")     # Minimum x-value of any turbine
        wt_all_x_max = -1*float("inf")
        wt_all_y_min = float("inf")
        wt_all_y_max = -1*float("inf")
        wt_all_z_max = -1*float("inf")  # Tallest rotor disk point of any turbine
        Drot_max = -1*float("inf")
        for turbkey in self.wts:
            wt_all_x_min = min(wt_all_x_min, self.wts[turbkey]['x'])
            wt_all_x_max = max(wt_all_x_max, self.wts[turbkey]['x'])
            wt_all_y_min = min(wt_all_y_min, self.wts[turbkey]['y'])
            wt_all_y_max = max(wt_all_x_min, self.wts[turbkey]['y'])
            wt_all_z_max = max(wt_all_z_max, self.wts[turbkey]['zhub'] + 0.5*self.wts[turbkey]['D'])
            Drot_max = max(Drot_max, self.wts[turbkey]['D'])
            
        xlow_lr_min  = wt_all_x_min - self.buffer_lr[0] * Drot_max
        xhigh_lr_max = wt_all_x_max + self.buffer_lr[1] * Drot_max 
        ylow_lr_min  = wt_all_y_min - self.buffer_lr[2] * Drot_max 
        yhigh_lr_max = wt_all_y_max + self.buffer_lr[3] * Drot_max 
        zhigh_lr_max = wt_all_z_max + self.buffer_lr[4] * Drot_max 

        # Calculate the minimum/maximum LR domain coordinate lengths & number of grid cells
        xdist_lr_min = xhigh_lr_max - xlow_lr_min  # Minumum possible length of x-extent of LR domain
        self.xdist_lr = self.ds_low_les * np.ceil(xdist_lr_min/self.ds_low_les)  # The `+ ds_lr` comes from the +1 to NS_LOW in Sec. 4.2.15.6.4.1.1
        # TODO: adjust xdist_lr calculation by also using `inflow_deg`
        self.nx_lr = int(self.xdist_lr/self.ds_low_les) + 1

        ydist_lr_min = yhigh_lr_max - ylow_lr_min
        self.ydist_lr = self.ds_low_les * np.ceil(ydist_lr_min/self.ds_low_les)
        # TODO: adjust ydist_lr calculation by also using `inflow_deg`
        self.ny_lr = int(self.ydist_lr/self.ds_low_les) + 1

        self.zdist_lr = self.ds_low_les * np.ceil(zhigh_lr_max/self.ds_low_les)
        self.nz_lr = int(self.zdist_lr/self.ds_low_les) + 1

        ## Calculate actual LR domain extent
        #   NOTE: Sampling planes should measure at AMR-Wind cell centers, not cell edges
        #   NOTE: Should we use dx/dy/dz values here or ds_lr?
        #           - AR: I think it's correct to use ds_lr to get to the xlow values,
        #               but then offset by 0.5*amr_dx0 if need be
        self.xlow_lr = self.ds_low_les * np.floor(xlow_lr_min/self.ds_low_les) - 0.5*self.dx0 + self.prob_lo[0]%self.ds_low_les
        self.xhigh_lr = self.xlow_lr + self.xdist_lr
        self.ylow_lr = self.ds_low_les * np.floor(ylow_lr_min/self.ds_low_les) - 0.5*self.dy0 + self.prob_lo[1]%self.ds_low_les
        self.yhigh_lr = self.ylow_lr + self.ydist_lr
        self.zlow_lr = 0.5 * self.dz0  # Lowest z point is half the height of the lowest grid cell
        self.zhigh_lr = self.zlow_lr + self.zdist_lr
        self.zoffsets_lr = np.arange(self.zlow_lr, self.zhigh_lr+self.ds_low_les, self.ds_low_les) - self.zlow_lr

        ## Check domain extents
        if self.xhigh_lr > self.prob_hi[0]:
            raise ValueError(f"LR domain point {self.xhigh_lr} extends beyond maximum AMR-Wind x-extent!")
        if self.xlow_lr < self.prob_lo[0]:
            raise ValueError(f"LR domain point {self.xlow_lr} extends beyond minimum AMR-Wind x-extent!")
        if self.yhigh_lr > self.prob_hi[1]:
            raise ValueError(f"LR domain point {self.yhigh_lr} extends beyond maximum AMR-Wind y-extent!")
        if self.ylow_lr < self.prob_lo[1]:
            raise ValueError(f"LR domain point {self.ylow_lr} extends beyond minimum AMR-Wind y-extent!")
        if self.zhigh_lr > self.prob_hi[2]:
            raise ValueError(f"LR domain point {self.zhigh_lr} extends beyond maximum AMR-Wind z-extent!")
        if self.zlow_lr < self.prob_lo[2]:
            raise ValueError(f"LR domain point {self.zlow_lr} extends beyond minimum AMR-Wind z-extent!")

        ## Check grid placement
        self._check_grid_placement()

        ## Save out info for FFCaseCreation
        self.extent_low = self.buffer_lr

    def _check_grid_placement(self):
        '''
        Check the values of parameters that were calculated by _calc_sampling_params
        '''

        ## Check that sampling grids are at cell centers
        # Low resolution grid
        amr_xgrid_level0 = np.arange(self.prob_lo[0], self.prob_hi[0], self.dx0)
        amr_ygrid_level0 = np.arange(self.prob_lo[1], self.prob_hi[1], self.dy0)
        amr_zgrid_level0 = np.arange(self.prob_lo[2], self.prob_hi[2], self.dz0)

        amr_xgrid_level0_cc = amr_xgrid_level0 + 0.5*self.dx0  # Cell-centered AMR-Wind x-grid
        amr_ygrid_level0_cc = amr_ygrid_level0 + 0.5*self.dy0
        amr_zgrid_level0_cc = amr_zgrid_level0 + 0.5*self.dz0
        
        sampling_xgrid_lr = self.xlow_lr + self.ds_lr*np.arange(self.nx_lr)
        sampling_ygrid_lr = self.ylow_lr + self.ds_lr*np.arange(self.ny_lr)
        sampling_zgrid_lr = self.zlow_lr + self.zoffsets_lr

        # TODO: These for loops could be replaced with a faster operation
        for coord in sampling_xgrid_lr:
            if coord not in amr_xgrid_level0_cc:
                raise ValueError("Low resolution x-sampling grid is not cell cenetered with AMR-Wind's grid!")
        for coord in sampling_ygrid_lr:
            if coord not in amr_ygrid_level0_cc:
                raise ValueError("Low resolution y-sampling grid is not cell cenetered with AMR-Wind's grid!")
        for coord in sampling_zgrid_lr:
            if coord not in amr_zgrid_level0_cc:
                raise ValueError("Low resolution z-sampling grid is not cell cenetered with AMR-Wind's grid!")

        # High resolution grids (span the entire domain to make this check easier)
        amr_xgrid_refine = np.arange(self.prob_lo[0], self.prob_hi[0], self.dx_refine)
        amr_ygrid_refine = np.arange(self.prob_lo[1], self.prob_hi[1], self.dy_refine)
        amr_zgrid_refine = np.arange(self.prob_lo[2], self.prob_hi[2], self.dz_refine)

        amr_xgrid_refine_cc = amr_xgrid_refine + 0.5*self.dx_refine
        amr_ygrid_refine_cc = amr_ygrid_refine + 0.5*self.dy_refine
        amr_zgrid_refine_cc = amr_zgrid_refine + 0.5*self.dz_refine

        for turbkey in self.hr_domains:
            nx_hr = self.hr_domains[turbkey]['nx_hr']
            ny_hr = self.hr_domains[turbkey]['ny_hr']
            xlow_hr = self.hr_domains[turbkey]['xlow_hr']
            ylow_hr = self.hr_domains[turbkey]['ylow_hr']

            sampling_xgrid_hr = xlow_hr + self.ds_hr*np.arange(nx_hr)
            sampling_ygrid_hr = ylow_hr + self.ds_hr*np.arange(ny_hr)
            sampling_zgrid_hr = self.hr_domains[turbkey]['zlow_hr'] + self.hr_domains[turbkey]['zoffsets_hr']

            # TODO: These for loops could be replaced with a faster operation
            for coord in sampling_xgrid_hr:
                if coord not in amr_xgrid_refine_cc:
                    raise ValueError("High resolution x-sampling grid is not cell cenetered with AMR-Wind's grid!")
            for coord in sampling_ygrid_hr:
                if coord not in amr_ygrid_refine_cc:
                    raise ValueError("High resolution y-sampling grid is not cell cenetered with AMR-Wind's grid!")
            for coord in sampling_zgrid_hr:
                if coord not in amr_zgrid_refine_cc:
                    raise ValueError("High resolution z-sampling grid is not cell cenetered with AMR-Wind's grid!")

    def write_sampling_params(self, outdir=None):
        '''
        Write out text that can be used for the sampling planes in an 
          AMR-Wind input file

        outdir: str
            Input file
        '''
        outfile = os.path.join(outdir, 'sampling_config.i')
        if not os.path.exists(outdir):
            print(f'Path {outdir} does not exist. Creating it')
            os.makedirs(outdir)
        if os.path.isfile(outfile):
            raise FileExistsError(f"{str(outfile)} already exists! Aborting...")

        print(f"Writing to {outfile} ...")
        with open(outfile,"w") as out:
            # Write high-level info for sampling
            sampling_labels_lr_str = " ".join(str(item) for item in self.sampling_labels_lr)
            sampling_labels_hr_str = " ".join(str(item) for item in self.sampling_labels_hr)
            out.write(f"# Sampling info generated by AMRWindSamplingCreation.py\n")
            out.write(f"incflo.post_processing                = {self.postproc_name_lr} {self.postproc_name_hr} # averaging\n\n")
            out.write(f"{self.postproc_name_lr}.output_format    = netcdf\n")
            out.write(f"{self.postproc_name_lr}.output_frequency = {self.output_frequency_lr}\n")
            out.write(f"{self.postproc_name_lr}.fields           = velocity # temperature tke\n")
            out.write(f"{self.postproc_name_lr}.labels           = {sampling_labels_lr_str}\n\n")

            out.write(f"{self.postproc_name_hr}.output_format    = netcdf\n")
            out.write(f"{self.postproc_name_hr}.output_frequency = {self.output_frequency_hr}\n")
            out.write(f"{self.postproc_name_hr}.fields           = velocity # temperature tke\n")
            out.write(f"{self.postproc_name_hr}.labels           = {sampling_labels_hr_str}\n")

            # Write out low resolution sampling plane info
            zoffsets_lr_str = " ".join(str(int(item)) for item in self.zoffsets_lr)

            out.write(f"\n# Low sampling grid spacing = {self.ds_lr} m\n")
            out.write(f"{self.postproc_name_lr}.Low.type         = PlaneSampler\n")
            out.write(f"{self.postproc_name_lr}.Low.num_points   = {self.nx_lr} {self.ny_lr}\n")
            out.write(f"{self.postproc_name_lr}.Low.origin       = {self.xlow_lr:.1f} {self.ylow_lr:.1f} {self.zlow_lr:.1f}\n")  # Round the float output
            out.write(f"{self.postproc_name_lr}.Low.axis1        = {self.xdist_lr:.1f} 0.0 0.0\n")  # Assume: axis1 oriented parallel to AMR-Wind x-axis
            out.write(f"{self.postproc_name_lr}.Low.axis2        = 0.0 {self.ydist_lr:.1f} 0.0\n")  # Assume: axis2 oriented parallel to AMR-Wind y-axis
            out.write(f"{self.postproc_name_lr}.Low.normal       = 0.0 0.0 1.0\n")
            out.write(f"{self.postproc_name_lr}.Low.offsets      = {zoffsets_lr_str}\n")

            # Write out high resolution sampling plane info
            for turbkey in self.hr_domains:
                wt_x = self.wts[turbkey]['x']
                wt_y = self.wts[turbkey]['y']
                wt_D = self.wts[turbkey]['D']
                if 'name' in self.wts[turbkey].keys():
                    wt_name = self.wts[turbkey]['name']
                else:
                    wt_name = f'T{turbkey}'
                sampling_name = f"High{wt_name}_inflow0deg"
                nx_hr = self.hr_domains[turbkey]['nx_hr']
                ny_hr = self.hr_domains[turbkey]['ny_hr']
                xlow_hr = self.hr_domains[turbkey]['xlow_hr']
                ylow_hr = self.hr_domains[turbkey]['ylow_hr']
                zlow_hr = self.hr_domains[turbkey]['zlow_hr']
                xdist_hr = self.hr_domains[turbkey]['xdist_hr']
                ydist_hr = self.hr_domains[turbkey]['ydist_hr']
                zoffsets_hr = self.hr_domains[turbkey]['zoffsets_hr']
                zoffsets_hr_str = " ".join(str(int(item)) for item in zoffsets_hr)

                out.write(f"\n# Turbine {wt_name} at (x,y) = ({wt_x}, {wt_y}), with D = {wt_D}, grid spacing = {self.ds_hr} m\n")
                out.write(f"{self.postproc_name_hr}.{sampling_name}.type         = PlaneSampler\n")
                out.write(f"{self.postproc_name_hr}.{sampling_name}.num_points   = {nx_hr} {ny_hr}\n")
                out.write(f"{self.postproc_name_hr}.{sampling_name}.origin       = {xlow_hr:.1f} {ylow_hr:.1f} {zlow_hr:.1f}\n")  # Round the float output
                out.write(f"{self.postproc_name_hr}.{sampling_name}.axis1        = {xdist_hr:.1f} 0.0 0.0\n")  # Assume: axis1 oriented parallel to AMR-Wind x-axis
                out.write(f"{self.postproc_name_hr}.{sampling_name}.axis2        = 0.0 {ydist_hr:.1f} 0.0\n")  # Assume: axis2 oriented parallel to AMR-Wind y-axis
                out.write(f"{self.postproc_name_hr}.{sampling_name}.normal       = 0.0 0.0 1.0\n")
                out.write(f"{self.postproc_name_hr}.{sampling_name}.offsets      = {zoffsets_hr_str}\n")
