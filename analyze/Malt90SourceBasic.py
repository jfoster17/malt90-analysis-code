
import malt_params as malt

class Malt90SourceBasic:
    """A class to store basic information about a Malt90Source.
    Can run on Draco, so limited imports required.
    Stores information about
    (1) When source was observed
    (2) What the noise properties are like
    (3) Which files went into the reduction. 
    """
    def __init__(self,atlasgal_id):
        self.id = atlasgal_id
        self.comments = self.get_comments()
        self.get_log_info()

    def get_log_info(self):
        self.location = malt.log_location
        f = open(self.location,'r')
        lines = f.readlines()
        self.raw_filenames = []
        self.source_names  = []
        for line in lines:
            if line.startswith("#"):
                pass
            else:
                en = line.split()
                try:
                    raw_filename = en[0]
                    atlasgal_id  = en[1]
                    sourcename   = en[2]
                    if atlasgal_id == self.id:
                        self.raw_filenames.append(raw_filename)
                        self.source_names.append(sourcename)
                except IndexError: #Deal with blank line
                    print("Blank or corrupt line in reduction log")
        f.close()

    def get_comments(self):
        """Get freeform comments about a source from text database."""
        pass
