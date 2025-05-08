import os
from time import sleep
from loky import get_reusable_executor


def say_hello(k):
    pid = os.getpid()
    print("Hello from {} with arg {}".format(pid, k))
    sleep(5.0)

    return k*k

executor = get_reusable_executor(max_workers=4, timeout=2)

input_list=[1,2,3,4,5,6,7,8]

results = list(executor.map(say_hello, input_list))
print(results)