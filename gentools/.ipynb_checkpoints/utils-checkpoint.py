import yaml

class Configfile:
    '''Class that hold configfile'''
    def __init__(self, configfile):
        self._yaml = self.__read_yaml(configfile)
    
    def __read_yaml(self, file):
        with open(file, 'r') as f:
            return yaml.full_load(f)
        
    def __getitem__(self, item):
        return self._yaml[item]
