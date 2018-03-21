#!/usr/bin/python3

import sys
import math as mt
import bisect
from functools import partial
import itertools
from operator import itemgetter
from itertools import groupby
def optparser(argv=sys.argv[1:]):
	
	usage = ''' 
	Usage:

	 python3 prepare_bias_1.6.py -m receptor.XX.map -b bias_file.dat -o receptor.XX.bias.map 

		-m : original map
		-b : tab separated input file (X Y Z Energy/PFP R90)
		-o : modified map
	'''

	opt = ('-m', '-b', '-o')

	if len(argv) == 6:
		args = { a:argv[i+1] for i, a in enumerate(argv) if a in opt }
	else:
		print('Error: invalid number of arguments')
		print(usage) 
		sys.exit(2)
	
	for i in opt:

		if i not in args:
			print(i, 'not found.')			
			print(usage)
			sys.exit(2) 
		
		if i == '-m' and args[i][-4:] != '.map':
			print(args['-m'], 'is not XX.map file.')
			sys.exit(2)
		
		if i == '-b' and args[i][-4:] != '.dat':
			print(args['-b'], 'is not bias dat file.')
			sys.exit(2)

	return {'map_in': args['-m'], 'map_out': args['-o'], 'bias_file': args['-b']}

class Grid():

	''' returns a dictionary with the parameters with which a grid was built
	 '''
		
	def __init__(self, map_in):
		
		self.grid = self.parse_file(map_in)
		self.grid = self.adicionales(self.grid)

	def parse_file(self, map_in):

		grid = dict()

		with open(map_in) as f:
			
			for _ in range(6):
				
				line = f.readline().split()
				
				if line[0] == 'SPACING':
					grid.update( {'spacing' : float(line[1]) })
					
				if line[0] == 'NELEMENTS':
					grid.update( {'NX' : int(line[1]) + 1 })
					grid.update( {'NY' : int(line[2]) + 1 })
					grid.update( {'NZ' : int(line[3]) + 1 })

				if line[0]=='CENTER':
					grid.update( {'CX' : float(line[1]) })
					grid.update( {'CY' : float(line[2]) })
					grid.update( {'CZ' : float(line[3]) })

				if line[0] == 'MACROMOLECULE':
					grid.update( {'file' : line[1] })
		return(grid)
	
	def adicionales(self, grid):

		grid.update( {'Mx' : grid['CX'] +  ((grid['NX'] - 1 ) / 2 ) * grid['spacing'] })
		grid.update( {'My' : grid['CY'] +  ((grid['NY'] - 1 ) / 2 ) * grid['spacing'] })
		grid.update( {'Mz' : grid['CZ'] +  ((grid['NZ'] - 1 ) / 2 ) * grid['spacing'] })
		
		grid.update( {'mx' : grid['CX'] -  ((grid['NX'] - 1 ) / 2 ) * grid['spacing'] })
		grid.update( {'my' : grid['CY'] -  ((grid['NY'] - 1 ) / 2 ) * grid['spacing'] })
		grid.update( {'mz' : grid['CZ'] -  ((grid['NZ'] - 1 ) / 2 ) * grid['spacing'] })
		
		self.grid.update( {'total_points': self.grid['NX'] * self.grid['NY'] * self.grid['NZ'] })

		return(grid)

class WaterSites():
	""" from bias file returns a dictionary with watersites properties
	"""
		
	def __init__(self, bias_file):

		self.watersites = self.parse_file(bias_file)
	
	def parse_file(self, bias_file):

		with open(bias_file) as f:
			
			WS = []

			for i, line in enumerate(f.readlines()):
				#if i != 0:
				line = line.split()
				
				if ((len(line) > 4) and ((line[0].split('.')[0]).strip("-").isnumeric())):
					ws = {	'x':float(line[0]),
							'y':float(line[1]),
							'z':float(line[2]),
							'R90':float(line[4].lstrip('\n')) }
					# Energy bias definition
					if float(line[3]) < 0:
						ws.update( {'dE':float(line[3]) } )
					# PFP bias definition
					else:
						ws.update( {'WFP':float(line[3])} )
								
					WS.append(ws)
				
		return(WS)

class Cubes():
	
	def __init__(self, watersites, grid):
		
		cubos = [self.obtener_cubo(ws, grid) for ws in watersites]

		self.cubos = [cubo for cubo in cubos if self.validar_cubo(cubo, grid)]

	def obtener_cubo(self, ws, grid):

		puntos_desde_centro = mt.ceil(ws['R90'] * 2 / grid['spacing'])

		costado_cubo = puntos_desde_centro * 2  + 1 

		total_puntos_cubo = costado_cubo ** 3		
		
		cnx = round( (ws['x'] - grid['mx'] ) / grid['spacing'] )
		cny = round( (ws['y'] - grid['my'] ) / grid['spacing'] )
		cnz = round( (ws['z'] - grid['mz'] ) / grid['spacing'] )

		center_cube_x = cnx * grid['spacing'] + grid['mx']
		center_cube_y = cny * grid['spacing'] + grid['my']
		center_cube_z = cnz * grid['spacing'] + grid['mz']
		
		center_point_map = 6 +  cnx + ( grid['NX']  ) * cny +  ( grid['NX']  ) * ( grid['NY'] )  * cnz 
		
		min_cube_x = center_cube_x - puntos_desde_centro * grid['spacing']
		min_cube_y = center_cube_y - puntos_desde_centro * grid['spacing']
		min_cube_z = center_cube_z - puntos_desde_centro * grid['spacing']

		min_nx = cnx - puntos_desde_centro
		min_ny = cny - puntos_desde_centro
		min_nz = cnz - puntos_desde_centro

		cut_x = min_cube_x + ( costado_cubo - 1 ) * grid['spacing']
		cut_y = min_cube_y + ( costado_cubo - 1 ) * grid['spacing']
		cut_z = min_cube_z + ( costado_cubo - 1 ) * grid['spacing']

		distancia_a_grilla =( (ws['x'] -center_cube_x )**2 + (ws['y'] - center_cube_y)**2 + (ws['z'] - center_cube_z)**2 ) ** (1/2) 

		prop_cubo = {'limite_grilla' : {'NX': grid['NX'], 'NY': grid['NY'], 'NZ': grid['NZ']},
					 'watersite_origen': {'x': ws['x'], 'y': ws['y'], 'z': ws['z']},
					 'centro_coordenadas': {'x': center_cube_x, 'y': center_cube_y, 'z': center_cube_z},
					 'minimo_coordenadas': {'x': min_cube_x, 'y': min_cube_y, 'z': min_cube_z},
					 'centro_desde_origen_puntos': {'x': cnx, 'y': cny, 'z': cnz},
					 'minimo_desde_origen_puntos': {'x': min_nx, 'y': min_ny, 'z': min_nz },
					 'limite': {'x': cut_x, 'y': cut_y, 'z':cut_z},
					 'center_point_map': int(center_point_map +1),
					 'spacing': grid['spacing'],
					 'puntos_totales': total_puntos_cubo,
					 'R90': ws['R90'],
					 'distancia_a_grilla': distancia_a_grilla
					 }

		if 'dE' in ws: prop_cubo.update( {'dE': ws['dE'] } )
		elif 'WFP' in ws: prop_cubo.update( {'WFP': ws['WFP']} )
		else: print('no se encuentra WFP ni dE en el ws ', prop_cubo['watersite_origen'])
		
		return prop_cubo

	def validar_cubo(self, cubo, grid):
		'''se observa si los cubos entran en el espacio de la grilla'''

		valido = True

		if   grid['Mx'] < cubo['limite']['x'] : valido = False
		elif grid['My'] < cubo['limite']['y'] : valido = False
		elif grid['Mz'] < cubo['limite']['z'] : valido = False

		elif grid['mx'] > cubo['minimo_coordenadas']['x'] : valido = False
		elif grid['my'] > cubo['minimo_coordenadas']['y'] : valido = False
		elif grid['mz'] > cubo['minimo_coordenadas']['z'] : valido = False

		if not valido:
			print('W: Cube with center coordinate', cubo['centro_coordenadas'], 'points out grid.' )
			return False

		return True

class Insert_bias():

	""" routine that inserts Gaussian modifications

	"""

	def __init__(self, cubes, map_in, map_out):

		modificaciones_por_cubo = []

		for cubo in cubes:
			modificaciones_por_cubo.append(self.modificacion_cubo(cubo))
		#dictionary or list of energy modifications (#point, energy)
		modifications_lst=self.merge_boxes(modificaciones_por_cubo)
		#modify the original map with the ws information.
		self.insert_bias(cubes, modifications_lst, map_in, map_out)
		

	def modificacion_cubo(self, cubo):

		modificacion = []
		
		RT = 0.593

		x = cubo['minimo_coordenadas']['x']
		y = cubo['minimo_coordenadas']['y']
		z = cubo['minimo_coordenadas']['z']

		nx = cubo['minimo_desde_origen_puntos']['x']
		ny = cubo['minimo_desde_origen_puntos']['y']
		nz = cubo['minimo_desde_origen_puntos']['z']

		for _ in range(cubo['puntos_totales']):

			if x > cubo['limite']['x']:
					x =  cubo['minimo_coordenadas']['x'] 
					y += cubo['spacing']
					
					nx = cubo['minimo_desde_origen_puntos']['x']
					ny += 1


			if y > cubo['limite']['y']:
				y =  cubo['minimo_coordenadas']['y']
				z += cubo['spacing']
				
				ny = cubo['minimo_desde_origen_puntos']['y']
				nz += 1

			distancia = ( (x - cubo['watersite_origen']['x'])**2 + (y - cubo['watersite_origen']['y'])**2 + (z - cubo['watersite_origen']['z'])**2 ) ** (1/2)
			dist_centro_cubo =  ( (x - cubo['centro_coordenadas']['x'])**2 + (y - cubo['centro_coordenadas']['y'])**2 + (z - cubo['centro_coordenadas']['z'])**2 ) ** (1/2)
			

			if 'WFP' in cubo:
				dE = -1 * RT * mt.log(cubo['WFP']) * mt.exp (-1 * (distancia ** 2) / (cubo['R90'] ** 2))
				dEo = -1 * RT * mt.log(cubo['WFP'])
	
			elif 'dE' in cubo:
				dE = cubo['dE'] * mt.exp (-1 * (distancia ** 2) / (cubo['R90'] ** 2))	

			if dE < -0.01 :
				posicion_OA_map = 6 +  nx + ( cubo['limite_grilla']['NX']  ) * ny +  ( cubo['limite_grilla']['NX'] ) * ( cubo['limite_grilla']['NY'] )  * nz
			
				modificacion.append((posicion_OA_map, dE))


			x += cubo['spacing']
			nx += 1
					
		return modificacion


	def merge_boxes(self, modificaciones):
		if len(modificaciones) > 0:
			flat_list = list(itertools.chain.from_iterable(modificaciones))
			#sort dic by index and energy.
			flat_list.sort(key=itemgetter(0,1))
			#save the lower energy for each tuple.
			lst_unique=[(key,) + tuple(elem for _, elem in group) for key, group in groupby(flat_list, lambda pair: pair[0])]
			#Generate a diccionary (index, energy)
			energy_dictionary = {}
			for i in lst_unique:
				energy_dictionary[i[0]] = i[1]
		#return dictionary dic(Number_point,energy)
		return energy_dictionary

	def insert_bias(self, cubos, modificaciones, map_in, map_out):
		Eori = {}
		with open(map_in) as f:
			with open(map_out, 'w') as g:
				for index_map_in, original_E in enumerate(f.readlines()):
					#find point in a dictionary.
					if index_map_in in modificaciones:
						# Define original energy 
						Eori[index_map_in] = float(original_E.rstrip('\n'))
						# Define new energy with ws modification.
						E = float(original_E.rstrip('\n')) + modificaciones[index_map_in]
						# write new energy value.
						g.write('{:.3f}'.format(E) + '\n')
				
					else:
						#write old energy value.	
						g.write(original_E)
						
		# output filename ( .map)
		dat_modif = 'grid_modif_' + map_in.split('.')[-2] +'.dat'
		# output variable with ws information.
		lineas = []
		with open(dat_modif, 'w') as p:
			p.write("\t".join(['x', 'y', 'z', 'point', 'D_min', 'E_ori', 'E_modif', 'E_delta', 'R90'])+'\n'	)
			for cubo in cubos:
				linea = ["%0.3f" % cubo['centro_coordenadas']['x'],
						 "%0.3f" % cubo['centro_coordenadas']['y'], 
						 "%0.3f" % cubo['centro_coordenadas']['z'],
						 str(cubo['center_point_map']), 
						 "%0.3f" % cubo['distancia_a_grilla'],
						 "%0.3f" % Eori[cubo['center_point_map']], 
						 "%0.3f" % (float(Eori[cubo['center_point_map']]) + float(modificaciones[cubo['center_point_map']])), 
						 "%0.3f" % float(modificaciones[cubo['center_point_map']]), 
						 "%0.3f" % cubo['R90']]
				lineas.append("\t".join(linea))
			for linea in lineas:
				p.write(linea + '\n')
#	Main
archivos = optparser()

grilla = Grid( archivos['map_in'] )

aguas = WaterSites( archivos['bias_file'] )

cubos = Cubes( aguas.watersites, grilla.grid )

modificaciones = Insert_bias( cubos.cubos, archivos['map_in'], archivos['map_out'] )
# 
