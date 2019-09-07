'''
This module defines the tao_plot_dict class, which keeps track of what
plots have been placed in tao along with what region they have been placed in.
Comes with methods for adding and removing plots to/from the dictionary.
'''
class tao_plot_dict(dict):
    '''
    Specialized dictionary for keeping track of plots and regions in tao
    '''
    def __init__(self, pipe, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self.pipe = pipe
    def place_template(self, template, region=None):
        '''
        Places the given template plot in the specified region
        If no region is specified, automatically chooses the next available region
        Returns the name of the region where the plot was placed
        '''
        plot_list_r = self.pipe.cmd_in('python plot_list r').splitlines()
        index_list = []
        region_list = []
        for i in range(len(plot_list_r)):
            ix, region_name = plot_list_r[i].split(';')[:2]
            if region:
                if region not in [region_name, '@R'+ix]:
                    continue
            elif ('@R'+ix in self.keys()) or (region in self.keys()):
                continue
            self['@R'+ix] = template
            self[region_name] = template
            self.pipe.cmd_in('place -no_buffer ' + '@R'+ix + ' ' + template)
            #self.pipe.cmd_in('set plot ' + '@R'+ix + ' visible = T')
            return region_name
        # Try to place in any region if region was specified
        if region:
            self.place_template(template, None)
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
        plot_list_r = self.pipe.cmd_in('python plot_list r').splitlines()
        for line in plot_list_r:
            if region in line.split(';'):
                ix_region = '@R' + line.split(';')[0] #e.g. @R3
                name_region = line.split(';')[1] #e.g. r11, top, etc
                # remove from self
                if ix_region in self:
                    self.pop(ix_region)
                if name_region in self:
                    self.pop(name_region)
                # unplace in tao
                self.pipe.cmd_in('place -no_buffer ' + ix_region + ' none')
