#!/usr/bin/python

import sys, getopt, re, textwrap

pi = 3.141592654

class ele_struct:
	def __init__(self):
		self.name = ''
		self.type = ''
		self.params = dict()
		self.printed = False

inputfile = sys.argv[1]
outputfile = sys.argv[2]

print 'Input file is ', inputfile
print 'Output file is ', outputfile

fin = open(inputfile,'r')
fou = open(outputfile,'w')
fleftovers = open('sad_to_bmad.leftovers','w')

latticeDefiners = ("drift", "bend", "quad", "sext", "oct", "mult", "sol", "cavi", "moni", "line", "beambeam", "apert", "mark")
regexDef = r"(.*?)=(\(.*?\))"
regexLine = r"\w"
regexParam = r"(.*?)=(\S*\w\.? (deg)?)"

param = [None]*100  #declare param as list with 100 elements
value = [None]*100  #declare param as list with 100 elements

def fouWrapWrite(line):
	MAXLINELENGTH = 120
	lines = textwrap.wrap(line,MAXLINELENGTH)
	tab = False
	for line in lines[:-1]:
		if tab:
			fou.write('         '+line+' &\n')
		else:
			fou.write(line+' &\n')
			tab = True
	if tab:
		fou.write('         '+lines[-1]+'\n')
	else:
		fou.write(lines[-1]+'\n')

ele = []

def makeBMAD(type, name, nparams, param, value):
	ele.append(ele_struct())
	ele[-1].name = name	
	known = True
	if type == 'drift':
		ele[-1].type = 'drift'
	elif type == 'bend':
		ele[-1].type = 'rbend'
	elif type == 'quad':
		ele[-1].type = 'quadrupole'
	elif type == 'sext':
		ele[-1].type = 'sextupole'
	elif type == 'oct':
		ele[-1].type = 'octupole'
	elif type == 'mult':
		ele[-1].type = 'sad_mult'
	elif type == 'sol':
		ele[-1].type = 'solenoid'
	elif type == 'cavi':
		ele[-1].type = 'rfcavity'
	elif type == 'moni':
		ele[-1].type = 'marker'
	elif type == 'mark':
		ele[-1].type = 'marker'
	elif type == 'beambeam':
		ele[-1].type = 'marker'
	elif type == 'apert':
		ele[-1].type = 'marker'
	else:
		fleftovers.write('unknown type: '+type+'\n')
		known = False

	if(known):
		for ix in range(0,nparams):
			parm = param[ix]
			val  = value[ix]
			if parm == 'l':
				ele[-1].params['l']= val
			elif parm == 'dx':
				ele[-1].params['x_offset']= val
			elif parm == 'dy':
				ele[-1].params['y_offset']= val
			elif parm == 'dz':
				ele[-1].params['z_offset']= val
			elif parm == 'radius':
				ele[-1].params['aperture']= val
			elif parm == 'angle':
				ele[-1].params['angle']= val
			elif parm == 'e1':
				ele[-1].params['e1']= val
			elif parm == 'e2':
				ele[-1].params['e2']= val
			elif parm == 'rotate':
				if type == 'bend':
					val, _, units = val.partition(" ")
					if units.lower() == 'deg':
						val = str(float(val)*2.0*pi/360.0)
					ele[-1].params['ref_tilt']= val
				else:
					ele[-1].params['tilt']= val
			elif parm == 'bz':
				ele[-1].params['bs_field']= val
			elif parm == 'k1':
				if type != 'mult':
					ele[-1].params['k1']= val
				else:
					ele[-1].params['a1']= val
			elif parm == 'k2':
				if type != 'mult':
					ele[-1].params['k2']= val
				else:
					ele[-1].params['a2']= val
			elif parm in ('k0','k3','k4','k5','k6','k7','k8','k9','k10','k11','k12','k13','k14','k15','k16','k17','k18','k19','k20','k21'):
				ele[-1].params['a'+parm[-1:]]= val
			elif parm in ('sk0','sk1','sk2','sk3','sk4','sk5','sk6','sk7','sk8','sk9','sk10','sk11','sk12','sk13','sk14','sk15','sk16','sk17','sk18','sk19','sk20','sk21'):
				ele[-1].params['b'+parm[-1:]]= val
			elif parm == 'f1':
				ele[-1].params['f1']= val
			elif parm == 'f2':
				ele[-1].params['f2']= val
			elif parm == 'fringe':
				if val == '0':
					ele[-1].params['fringe_at'] = 'both_ends$'
					ele[-1].params['fringe_kind'] = 'nonlin_only$'
				elif val == '1':
					ele[-1].params['fringe_at'] = 'entrance_end$'
					ele[-1].params['fringe_kind'] = 'full$'
				elif val == '2':
					ele[-1].params['fringe_at'] = 'exit_end$'
					ele[-1].params['fringe_kind'] = 'full$'
				elif val == '3':
					ele[-1].params['fringe_at'] = 'both_ends$'
					ele[-1].params['fringe_kind'] = 'full$'
			elif parm == 'eps':
				ele[-1].params['eps_step_scale']= val
			elif parm == 'freq':
				ele[-1].params['rf_frequency']= val
			elif parm == 'volt':
				ele[-1].params['voltage']= val
			elif parm == 'bound':
				ele[-1].params['bound']= val
			elif parm == 'geo':
				ele[-1].params['geo']= val
			else:
				# chi1, chi2, chi3
				leftoverline = name
				leftoverline += ' >>>'+type+'<<< '
				leftoverline += ' >>>'+parm+'<<< '
				leftoverline += ' >>>'+val+'<<<\n'
				fleftovers.write(leftoverline)

lattice = []
def makeLINE(type,definitions):
	linename, _, latline = definitions.partition("=")
	lattice.append(linename)
	latline = latline.strip()
	latline = latline.lower()
	latline = latline[1:]
	latline = latline[:-1]
	for elename in re.split('\W+',latline)[:-1]:
		lattice.append(elename)

def parseDirective(directive):
	directive = directive.strip()  #Remove leading and trailing blanks.
	type, _, definitions = directive.partition(" ")
	if type == 'line':
		makeLINE(type,definitions)
	else:
		if type in latticeDefiners:
			for name_and_params in re.finditer(regexDef, definitions):
				name = name_and_params.group(1).strip()   #name is the name of the element
				params = name_and_params.group(2).strip() #params is a continuous string containing all the parameters for that element
				params = params[1:] 
				params = params[:-1]
				nparams = 0
				for param_and_val in re.finditer(regexParam, params):
					param[nparams] = param_and_val.group(1).strip()
					value[nparams] = param_and_val.group(2).strip()
					nparams += 1
				makeBMAD(type, name, nparams, param, value)

#Read in SAD file line-by-line.  Assemble lines into directives, which are delimited by a ; (colon).
#Call parseDirective whenever an entire directive has been obtained.
directive = ''
for line in fin:
	line = line.strip()    #Remove leading and trailing blanks.
	line = line.lower()    #All letters to lower case.
	if line:        #Detect blank line.
		if line[0] != '!':     #Detect comment line.
			if not ((line[0:2] == '(*') and (line[-2:] == '*)')):
				directive = directive + line + " "
				if directive[-2:] == '; ':  #End of directive 
					directive = directive[:-2]
					parseDirective(directive)
					directive = ''

inside_sol = False
for element in lattice[1:]:
	ix = [x.name for x in ele].index(element)
	if ele[ix].type == 'solenoid':
		bs_field = ele[ix].params['bs_field']	
		if 'bound' in ele[ix].params:
			if ele[ix].params['bound'] == '1':
				if inside_sol == False:
					inside_sol = True
				else:
					inside_sol = False
	else:
		if ele[ix].printed == False:
			bmadline = ele[ix].name+': '+ele[ix].type
			for param in iter(ele[ix].params):
				bmadline += ', '+param+'='+ele[ix].params[param]
			if inside_sol:
				bmadline += ', '+'bs_field='+bs_field
			fouWrapWrite(bmadline)
			ele[ix].printed = True

bmadline = lattice[0]+': line = ('
for name in lattice[1:]:
	bmadline += name+', '
bmadline = bmadline[:-2]+')'
fouWrapWrite(bmadline)

fin.close()
fou.close()
fleftovers.close()











#
