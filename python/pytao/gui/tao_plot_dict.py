'''
This module defines the tao_plot_dict class, which keeps track of what
plots have been placed in tao along with what region they have been placed in.
Comes with methods for adding and removing plots to/from the dictionary.
'''
class tao_plot_dict(dict):
    '''
    Specialized dictionary for keeping track of plots and regions in tao
    mode should be either "matplotlib" or "pgplot"
    '''
    def __init__(self, mode, pipe, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self.pipe = pipe
        self.mode = mode
        # Take note of the allowed region names and indices
        self._known_region_names = self.pipe.cmd_in("python plot_list r").splitlines()
        self._known_region_nums = [0]*len(self._known_region_names)
        for i in range(len(self._known_region_names)):
            self._known_region_names[i] = self._known_region_names[i].split(';')
            self._known_region_nums[i] = '@R' + self._known_region_names[i][0]
            self._known_region_names[i] = self._known_region_names[i][1]
        self.pgplot_settings = {"rows":"1", "columns":"1", "lat_layouts":"0"}

    def place_template(self, template, region=None):
        '''
        Places the given template plot in the specified region
        If no region is specified, automatically chooses the next available region
        Returns the name of the region where the plot was placed
        In pgplot mode, it is the caller's responsibility to ensure
        that a valid region name is provided, otherwise the plot will not be placed
        '''
        if self.mode == "matplotlib":
            return self._mpl_place_template(template, region)
        elif self.mode == "pgplot":
            return self._pgplot_place_template(template, region)

    def _pgplot_place_template(self, template, region):
        '''
        Implementation of place_template for pgplot mode
        Will only place the plot if a valid region name is provided
        '''
        if region:
            # Check that region is a valid region name
            if region in self._known_region_names:
                name = region
                num = self._known_region_names.index(region)
            elif region in self._known_region_nums:
                num = region
                name = self._known_region_nums.index(region)
            else: # region is not a valid region name
                return
            # Place the template
            self.pipe.cmd_in('place ' + name + ' ' + template)
            # Tao will unplace plots as necessary, so the
            # method below will update this dictionary's entries
            self.pgplot_update()

    def _mpl_place_template(self, template, region=None):
        '''
        Implementation of place_template specific to matplotlib mode
        '''
        # Get list of regions from tao
        plot_list_r = self.pipe.cmd_in('python plot_list r').splitlines()
        index_list = []
        region_list = []
        for i in range(len(plot_list_r)):
            ix, region_name = plot_list_r[i].split(';')[:2]
            if region:
                # Try to match region to region_name or @Rix
                if region not in [region_name, '@R'+ix]:
                    continue
            # Place in next available region otherwise
            elif ('@R'+ix in self.keys()) or (region in self.keys()):
                continue
            self['@R'+ix] = template
            self[region_name] = template
            self.pipe.cmd_in('place -no_buffer ' + '@R'+ix + ' ' + template)
            return region_name
        # Try to place in any region if region was specified
        if region:
            return self.place_template(template, None)


    def unplace_template(self, template):
        '''
        Remove the given template from self and unplace it in tao
        '''
        while template in self.values():
            for key, value in self.items():
                if value==template:
                    self.pop(key)
                    if key.find('@R')==0:
                        self.pipe.cmd_in('place -no_buffer ' + key + ' none')
                    break


    def unplace_region(self, region):
        '''
        Unplace whatever template is in the specified region, and remove the
        corresponding entry(ies) from self
        '''
        # Make sure region is a valid region name
        if region in self._known_region_names:
            name = region
            ix = self._known_region_names.index(region)
            num = self._known_region_nums[ix]
        elif region in self._known_region_nums:
            num = region
            ix = self._known_region_nums.index(region)
            name = self._known_region_names[ix]
        else:
            return # do nothing if region is not valid
        # Remove the plot at the specified region
        self.pipe.cmd_in('place -no_buffer ' + num + ' none')
        if name in self:
            self.pop(name)
        if num in self:
            self.pop(num)

    def pgplot_update(self):
        '''
        Reads the current state of plot regions in Tao, and
        sets self.pgplot_settings appropriately
        '''
        # Set the contents of this dictionary to match
        # the occupied pgplot regions
        self.clear()
        # Find the occupied regions
        plot_list_r = self.pipe.cmd_in("python plot_list r").splitlines()
        regions = []
        for line in plot_list_r:
            line = line.split(';')
            if line[3] == 'T':
                self['@R' + line[0]] = line[2]
                self[line[1]] = line[2]
                regions.append(line[1])
        # Determine the number of rows, columns, and lat_layouts
        r = 1
        c = 1
        lat = 0
        for region in regions:
            rcl = region_to_max_rcl(region)
            if rcl[0] > r:
                r = rcl[0]
            if rcl[1] > c:
                c = rcl[1]
            if rcl[2] > lat:
                lat = rcl[2]
        # Update self.pgplot_settings
        self.pgplot_settings["rows"] = str(r)
        self.pgplot_settings["columns"] = str(c)
        self.pgplot_settings["lat_layouts"] = str(lat)



#-----------------------------------------------------------
# Utilities
def region_to_max_rcl(region_name):
    '''
    Converts the region name (e.g. r2312) to a
    (row number, column number, lat_layout) tuple
    where the counts in the tuple are the maximum number of
    rows/columns/lat_layouts that are needed to accomadate
    the given region
    returns (-1,-1,-1) if the region_name is invalid
    '''
    # Case 1: top/bottom
    if region_name in ["top", "bottom"]:
        return (2,1,0)
    #Case 2: layout(1,2)
    if region_name[:6] == "layout":
        if len(region_name) == 6:
            lat = 1
        else:
            try:
                lat = int(region_name[6])
                lat = 2
            except ValueError:
                lat = -1
        if lat not in [1, 2]:
            return (-1,-1,-1)
        else:
            return (0,0,lat)
    # Case 3: rAB
    if len(region_name) == 3:
        try:
            r = int(region_name[1])
            rmax = int(region_name[2])
        except ValueError:
            return (-1, -1, -1)
        if (r not in range(1, rmax+1)) or (rmax<1) or (rmax>4):
            return (-1,-1,-1)
        else:
            if rmax==2:
                return (rmax,1,0)
            else:
                return (rmax,1,1)
    # Case 4: rABCD
    if len(region_name) == 5:
        try:
            r = int(region_name[3])
            rmax = int(region_name[4])
            c = 3 - int(region_name[1])
            cmax = int(region_name[2])
        except ValueError:
            return (-1, -1, -1)
        if ((r not in range(1, rmax+1)) or (rmax not in range(1,5))
                or (c not in range(1, cmax+1)) or (cmax not in range(1,3))):
            return (-1,-1,-1)
        else:
            return (rmax,cmax,1)

