
# coding: utf-8

# In[2]:


class Drug:
	def __init__(self,dbid,name,groups,classification,patent_date):
		self.dbid = dbid
		self.name = name
		self.edges = []
		self.classification = classification
		self.groups = groups
		self.patent_date = patent_date


class Protein:
	def __init__(self,upid,name,location,processes):
		self.upid = upid
		#self.organism = organism
		self.name = name
		self.subcell_location = location
		self.bio_processes = processes

