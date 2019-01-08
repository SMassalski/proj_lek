
# coding: utf-8

# In[2]:


class Drug:
	def __init__(self,dbid,name,groups,categories):
		self.dbid = dbid
		self.name = name
		self.edges = []
		self.categories = set(categories)
		self.groups = groups
		
	def add_edge(self,edge):
		self.edges.append(edge)
		


# In[3]:


class Edge:
	def __init__(self,drug,protein,subtype):
		self.drug = drug
		self.protein = protein
		self.subtype = subtype


# In[4]:


class Protein:
	def __init__(self,upid,organism):
		self.upid = upid
		self.organism = organism
		self.name = ''
		self.subcell_location = None
		self.bio_process = None

	def update(self,proc,loc):

		self.bio_process = set(proc.split('; '))
		self.subcell_location = set(loc.split('; '))

