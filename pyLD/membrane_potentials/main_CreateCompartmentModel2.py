#
#
#
import os, sys, glob
import numpy as np
import itertools
from pathlib import Path
import h5py
import pprint
import pickle
import trimesh

import networkx as nx
import matplotlib.pyplot as plt


class CreateCompartmentModel():

	def _loader(self):
		data_dir = r'C:\Users\uraku\Desktop\LatticeMicrobes\sim_membranepot\data'
		
		filename = os.path.join(data_dir, 'data_shared.pickle')
		with open(filename, 'rb') as f:
			self.org_graph = pickle.load(f)
			self.vertices  = pickle.load(f)

		filename = os.path.join(data_dir, 'data_spines.pickle')
		with open(filename, 'rb') as f:
			self.spines   = pickle.load(f)
			self.head_fs  = pickle.load(f)
			self.neck_fs  = pickle.load(f)

		filename = os.path.join(data_dir, 'data_dendrite.pickle')
		with open(filename, 'rb') as f:
			self.dendritic_nodes = pickle.load(f)
			self.dendritic_edges = pickle.load(f)

		self.nodes           = [n['node'] for n in self.dendritic_nodes]
		self.locations       = obtain_nodes_location(self.org_graph, self.nodes)
		self.container_faces = [n['faces'] for n in self.dendritic_edges]

	def obtain_id_node_from_old_node(self, id_old_node):
		ids = [ id for id, node in self.new_graph.nodes(data=True) if ('node_org' in node) and (node['node_org'] == id_old_node) ]
		return ids[0]

	def add_node(self, **kwargs):
		id = self.count_node
		self.new_graph.add_node(id, **kwargs)
		self.count_node += 1
		return id

	def add_edge(self, *args, **kwargs):
		id = self.count_edge
		kwargs.update({'id_compartment': id})
		self.new_graph.add_edge(*args, **kwargs)
		self.count_edge += 1
		return

	def add_neck(self, new_node_basement, neck):
		new_node_neck_head = self.add_node( attr = 'neck' )
		self.add_edge( new_node_neck_head, new_node_basement, \
			length   = neck['length'],   \
			volume   = neck['volume'],   \
			area     = neck['area'],     \
			id_neck  = neck['id_neck'],  \
			attr     = 'neck')
		return new_node_neck_head

	def add_neck_branched(self, new_node_basement, neck, j):
		new_node_neck = self.add_node( attr = 'neck' )
		self.add_edge( new_node_neck, new_node_basement, \
			length   = neck['lengths_for_branch'][j], \
			volume   = neck['volumes_for_branch'][j], \
			area     = neck['areas_for_branch'][j]  , \
			faces    = neck['faces_for_branch'][j]  , \
			attr     = 'neck-base')
		return new_node_neck

	def add_head(self, new_node_neck_head, head):
		# head
		node_org = head['node_terminal']
		loc      = self.org_graph.nodes[ node_org ]['loc']
		new_node_head = self.add_node( \
			node_org = node_org, \
			loc = loc, \
			attr = 'head-term' )
		self.add_edge( new_node_head, new_node_neck_head, \
			length   = head['length'],   \
			volume   = head['volume'],   \
			area     = head['area'],     \
			id_head  = head['id_head'],  \
			attr     = 'head'
			)

	def create(self):
		# Create a new graph.
		self.new_graph = nx.Graph()
		self.count_node = 0
		self.count_edge = 1 # for biochemical simulation

		# Add dendrite.
		new_nodes_dendrite = []
		for node in self.dendritic_nodes:
			node = self.add_node(\
				loc          = node['loc'], \
				node_org     = node['node'], \
				attr         = 'dendrite')
			new_nodes_dendrite.append(id)

		for edge in self.dendritic_edges:
			id_node0 = self.obtain_id_node_from_old_node( edge['edge'][0] )
			id_node1 = self.obtain_id_node_from_old_node( edge['edge'][1] )
			
			self.add_edge(id_node0, id_node1, \
				edge_org = edge['edge']				, \
				length   = edge['length']			, \
				volume   = edge['volume']			, \
				area     = edge['area']				, \
				faces    = edge['faces']			, \
				attr     = 'dendrite')
			# print('self.count_edge  ', self.count_edge)


		# Add spines.
		for spine in  self.spines :
			new_node_dend_neck = self.obtain_id_node_from_old_node( spine['node_dend'] )
			for neck, head in zip( spine['neck'], spine['head'] ):
				if 'area' in neck:
					# neck
					new_node_neck_head = self.add_neck(new_node_dend_neck, neck)
					# head
					self.add_head(new_node_neck_head, head)

				elif 'areas_for_branch' in neck:
					# neck(base)
					new_node_neck_base = self.add_neck_branched(new_node_dend_neck, neck, -1)
						
					for j in range( len(neck['areas_for_branch']) - 1 ):
						# neck 
						new_node_neck_head = self.add_neck_branched(new_node_neck_base, neck, j)
						# head
						self.add_head(new_node_neck_head, head[j])

	def __init__(self):
		self._loader()
		self.create()


if __name__ == '__main__':
	model = CreateCompartmentModel()


	filename = 'compartment_model.pickle'
	with open(filename, 'wb') as f:
		pickle.dump(model.new_graph, f)
		pickle.dump(model.vertices, f)
		pickle.dump(model.head_fs, f)
		pickle.dump(model.neck_fs, f)


