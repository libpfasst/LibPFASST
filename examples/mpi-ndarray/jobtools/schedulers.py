
from hopper import HopperPBS
from juqueen import JQLL
from serial import Serial

by_name = {
    'hopper': HopperPBS,
    'edison': HopperPBS,
    'juqueen': JQLL,
    'serial': Serial,
}
