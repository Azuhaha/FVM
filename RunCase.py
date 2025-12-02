from RunFalseTransientCase import *
from RunTransientCase import *
def runcase(reg):
    if reg.STEADY_STATE_RUN:
        runFalseTransientCase(reg)
    else:
        runTransientCase(reg)