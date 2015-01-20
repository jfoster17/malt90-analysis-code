import Malt90SourceBasic

class Malt90SourceSNR(Malt90SourceBasic.Malt90SourceBasic):

    def __init__(self,name):
        self.basic_name=name
        self.name = name
        self.data_dir = "/DATA/MALT_1/MALT90/data/"
        self.data_dirs = [self.data_dir+"year1/sources/"+self.basic_name+"/",self.data_dir+"year2/sources/"+self.basic_name+"/",self.data_dir+"year3/sources/"+self.basic_name+"/"]
