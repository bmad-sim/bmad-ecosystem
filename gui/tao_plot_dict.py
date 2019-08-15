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
      if region and (region not in [region_name, '@R'+ix]):
          continue
      elif ('@R'+ix in self.keys()) or (region in self.keys()):
        continue
      self['@R'+ix] = template
      self[region_name] = template
      self.pipe.cmd_in('place -no_buffer ' + '@R'+ix + ' ' + template)
      self.pipe.cmd_in('set plot ' + '@R'+ix + ' visible = T')
      return region_name
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
