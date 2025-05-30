"""Primary class for converting experiment-specific behavior."""
from pynwb.file import NWBFile

from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict

class Embargo2002ABehaviorInterface(BaseDataInterface):
    """Behavior interface for embargo2002a conversion"""

    keywords = ["behavior"]
    
    def __init__(self):
        # This should load the data lazily and prepare variables you need
        pass

    def get_metadata(self) -> DeepDict:
        # Automatically retrieve as much metadata as possible from the source files available
        metadata = super().get_metadata()   

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: dict):
        # All the custom code to add the data the nwbfile

        raise NotImplementedError()
