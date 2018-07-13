import inspect, os, sys
currentdir = os.path.dirname(os.path.abspath(os.path.abspath(inspect.getfile(inspect.currentframe()))))
for i in range(1):
    currentdir = os.path.dirname(currentdir)

if currentdir not in sys.path:
    sys.path.insert(0, currentdir)
