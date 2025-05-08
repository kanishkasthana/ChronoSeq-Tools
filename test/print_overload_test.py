from __future__ import print_function
import builtins as __builtin__
from datetime import datetime,timedelta

logFile=open('correctSubstitutionErrors.log', 'w')
def print(*args, **kwargs):
    output_string=""
    for arg in args:
        output_string+=str(arg)+" "
    logFile.write(output_string)
    return __builtin__.print(*args, **kwargs)

print("test",10,"what are you doing",datetime.now())