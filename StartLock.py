#Author Kanishk Asthana kasthana@eng.ucsd.edu
from multiprocessing.managers import BaseManager
from threading import Lock
import socket
lock=Lock()
print("Starting Global Read/Write Lock on :",socket.gethostname())

class LockManager(BaseManager):
    pass

LockManager.register("getLock", callable=lambda:lock)

manager = LockManager(address=(socket.gethostname(), 55442), authkey=b'Lock')

server = manager.get_server()
print("IP,Port: ",server.address)
server.serve_forever()