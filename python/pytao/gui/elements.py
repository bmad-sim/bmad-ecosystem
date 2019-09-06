from parameters import tao_parameter
from parameters import tao_parameter_dict

class lat_element():
    '''
    Holds the essential information for a given lattice element.
    A lat_element has the following properties:
    self.id: the ele identifier in the form uni@branch>>ele_ix|which
    self.params: An ordered dictionary of tao_parameters
        This holds the general parameters related to the element
    self.has: An ordered dictionary of bools
        This holds whether or not this element has a given property
    self.num: An ordered dictionary of ints
        This holds the number of something this element has
    '''
    def __init__(self, u_ix, ix_branch, ix_ele, which, pipe):
        # u_ix: universe index
        # ix_branch: branch index
        # ix_ele: element index
        # which: model, base, or design
        # pipe: the tao_interface object to use
        # ALL ARGUMENTS (besides pipe) MUST BE STRINGS
        self.pipe = pipe
        self.id = str(u_ix) + '@' + str(ix_branch) \
                + '>>' + str(ix_ele) + '|' + which
        fetch_cmd = "python ele:head " + self.id
        data_list = self.pipe.cmd_in(fetch_cmd)
        data_list = data_list.splitlines()

        p_list = [] # For parameters
        h_list = [] # For has#...
        n_list = [] # For num#...
        # All elements have gen_attribs:
        h_list.append("gen_attribs;LOGIC;F;T")
        for item in data_list:
            if item.find("has#") != -1:
                h_list.append(item[4:])
            elif item.find("num#") != -1:
                n_list.append(item[4:])
            else:
                # Manually set can_vary to F for base and design
                if which != "model":
                    item_parts = item.split(';')
                    item = ""
                    for i in range(len(item_parts)):
                        if i != 2:
                            item += item_parts[i] + ';'
                        else:
                            item += 'F;'
                    # Remove extra ;
                    item = item[:-1]
                p_list.append(item)

        self.params = tao_parameter_dict(p_list)
        self.has = tao_parameter_dict(h_list)
        self.num = tao_parameter_dict(n_list)
