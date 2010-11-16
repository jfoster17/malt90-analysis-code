
import glob
from asap import *
import os
import fcntl
import string
from subprocess import *

class ReduceLog:
	def __init__(self):
		self.location = '/DATA/MALT_1/MALT90/reduction_log.txt'
		self.lock = '/DATA/MALT_1/MALT90/lock.txt'
#		self.srcdir = '/DATA/ARCHIVE_1/Mopra/'		
		self.sd = '/nfs/atapplic/malt/reduce/'
		self.source_dir = '/DATA/MALT_1/MALT90/raw_data/'
		self.class_file = self.sd+'Malt90Catalog_classifications_v7.txt'
		self.fname = []
		self.id = []
		self.source = []
		self.rename = []
		self.smooth = []
		self.ldata = []
		self.gzilla = []
		self.arrange = []
		self.mommaps = []
	
	def set_val(self,source_file,field,value):
		g = open(self.lock)
		fcntl.flock(g,2) #Lock the lock file during access
		self.read()
		try:
			i = self.source.index(source_file)
		except ValueError:
			i = -1 # no match
		if field == "rename":
			self.rename[i] = value
		elif field == "smooth":
			self.smooth[i] = value
		elif field == "ldata":
			self.ldata[i] = value
		elif field == "gzilla":
			self.gzilla[i] = value
		elif field == "arrange":
			self.arrange[i] = value
		elif field == "mommaps":
			self.mommaps[i] = value
		self.save()
		os.system('''awk '{print $1",",$2",",$3",",$4",",$5",",$6",",$7",",$8",",$9}' /DATA/MALT_1/MALT90/reduction_log.txt > /DATA/MALT_1/MALT90/input_reduction_log.txt''')
		fcntl.flock(g,8) #Release lock file

	def check_val(self,source_file,field,vcheck):
		g = open(self.lock)
		fcntl.flock(g,2) #Lock the lock file during access
		self.read()
		try:
			i = self.source.index(source_file)
    		except ValueError:
        		i = -1 # no match
		if field == "rename":
			value = self.rename[i]
		elif field == "smooth":
			value = self.smooth[i]
		elif field == "ldata":
			value = self.ldata[i]
		elif field == "gzilla":
			value = self.gzilla[i]
		elif field == "arrange":
			value = self.arrange[i]
		elif field == "mommaps":
			value = self.mommaps[i]
		if float(value) < float(vcheck):
			return(True)
		else:
			return(False)
		fcntl.flock(g,8) #Release lock file

	def find_new_files(self,vcheck,datestring="20"):
		g = open(self.lock)
		fcntl.flock(g,2) #Lock the lock file during access
		self.read()
		new_files = []
#		print(self.rename)
		for i,entry in enumerate(self.rename):
		#	print(entry)
			
			if (entry < vcheck or self.smooth[i] < vcheck) and float(self.id[i]) > 0: #The ID check removes cal sources
				#print(self.id[i])
				#print(self.fname[i])
				if datestring in self.fname[i]:
					new_files.append(self.fname[i])
		fcntl.flock(g,8) #Release lock file
		return(new_files)

	def find_files_with_date(self,date):
		g = open(self.lock)
		fcntl.flock(g,2) #Lock the lock file during access
		self.read()
		new_files = []
		for i,entry in enumerate(self.fname):
			if date in entry and float(self.id[i]) > 0 and not 'B' in self.source[i]:
				ss = self.source[i].rstrip('GLatLon_')
				if ss not in new_files:
					new_files.append(ss)
		fcntl.flock(g,8) #Release lock file
		return(new_files)

	def get_name(self,filename):
		g = open(self.lock)
		fcntl.flock(g,2)
		self.read()
		i = self.fname.index(filename)
		fcntl.flock(g,8)
		return(self.id[i],self.source[i])		
	
	def update(self,date="all",update=True,in_middle_of_obs=False):
		if date == "all":
			raw_files = glob.glob(self.source_dir+'*M516.rpf')
		else:
			raw_files = glob.glob(self.source_dir+date+'*M516.rpf')
		raw_files.sort()
		g = open(self.lock)
		fcntl.flock(g,2) #Lock the lock file during access
		if update:
			self.read()
		for file in raw_files:
			if os.path.basename(file) in self.fname:
				print("Skipping file..."+os.path.basename(file))
			else:
				s = scantable(file,average=False)	
				possible_names = s.get_sourcename()
				try:
					source_name = possible_names[0].rstrip('_R')
				except IndexError:
					source_name = '--'
				source = ''
				atlasgal_id = '0'
				scan_dir = ''
				fullname_catalog =''
				try:
					(atlasgal_id,scan_dir) = source_name.split('_')
					atlasgal_id = atlasgal_id.lstrip('0')
				except:
					if source_name == 'g301cal':
						source = 'g301cal_____________'
					elif source_name == 'g337cal':
						source = 'g337cal_____________'
					else:
						print("Failed to find source name for "+file)
						source = '____________________'
				#Special case for badly-named 7mm pointing sources
				if atlasgal_id in string.letters:
					atlasgal_id = 0
					source = '____________________'

				f = open(self.class_file,'r')
				for line in f:
					pieces = line.split()
					atlasgal_catalog = pieces[0]
					if atlasgal_id == atlasgal_catalog:
						fullname_catalog = pieces[1]
						break
				f.close()
				st = os.stat(file)
				size = st[6] #Grab filesize in bytes
				if source:
					if size < 80000000:
						source = source+"B"
				if not source: #Not a calibration source. i.e. a real map
					if size > 300000000:
						#Give each valid entry a unique name
						source = fullname_catalog+'_'+scan_dir
						i = 2
						while source in self.source:
							source = fullname_catalog+'_'+str(i)+'_'+scan_dir
							i = i+1
					else:
						source = fullname_catalog+'_'+scan_dir+'B'	
					
				if source.endswith('B') and in_middle_of_obs: #Don't enter a bad source if we're in the middle of observing
					pass
					
				else:
					self.fname.append(os.path.basename(file))
					self.id.append(atlasgal_id)
					self.source.append(source)
					self.rename.append("0")
					self.smooth.append("0")
					self.ldata.append("0")
					self.gzilla.append("0")
					self.arrange.append("0")
					self.mommaps.append("0")
				print(os.path.basename(file)+'\t'+source)
		self.save()			
		os.system('''awk '{print $1",",$2",",$3",",$4",",$5",",$6",",$7",",$8",",$9}' /DATA/MALT_1/MALT90/reduction_log.txt > /DATA/MALT_1/MALT90/input_reduction_log.txt''')
		fcntl.flock(g,8) #Release lock file

	def read(self):
		f = open(self.location,'r')
		lines = f.readlines()
		self.fname = []
		self.id = []
		self.source = []
		self.rename = []
		self.smooth = []
		self.ldata = []
		self.gzilla = []
		self.arrange = []
		self.mommaps = []
		for line in lines:
			line.strip()
			if line.startswith("#"):
				pass
			else:
				en = line.split()
				self.fname.append(en[0])
				self.id.append(en[1])
				self.source.append(en[2])
				self.rename.append(en[3])
				self.smooth.append(en[4])
				self.ldata.append(en[5])
				self.gzilla.append(en[6])
				self.arrange.append(en[7])
				self.mommaps.append(en[8])
		f.close()

	def save(self):
		header = """#File			        ID      Name  	 	                Rename  Smooth  Ldata  Gzilla  Arrange  MomMaps
"""
		f = open(self.location,'w')
		f.write(header)
		decorated = zip(self.fname,self.id,self.source,self.rename,self.smooth,self.ldata,self.gzilla,self.arrange,self.mommaps)
		decorated.sort()
		self.fname,self.id,self.source,self.rename,self.smooth,self.ldata,self.gzilla,self.arrange,self.mommaps = zip(*decorated)
		for i,junk in enumerate(self.fname):
			outstring = self.fname[i]+'\t'+self.id[i]+'\t'\
			+self.source[i]+'\t\t'+self.rename[i]+'\t'+self.smooth[i]\
			+'\t'+self.ldata[i]+'\t'+self.gzilla[i]+'\t'\
			+self.arrange[i]+'\t'+self.mommaps[i]+'\n'
			f.write(outstring)		