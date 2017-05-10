from galaxy.datatypes import data
from galaxy.datatypes.binary import Binary

class Cel( Binary ):
    """Class for generic CEL binary format"""
    file_ext = "cel"
Binary.register_unsniffable_binary_ext("cel")