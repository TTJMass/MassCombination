# Runs doFit.py unchanged (blinded), records the x/y of plt.plot calls (chi2 scan points, blinded x) to scan_capture.json.
import sys, json, runpy, atexit, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = '/Users/sewuchte/Claude_MassComb/MassCombination/mtpole-ttj-pyconvino'
REC = []; _p = plt.plot
def p(*a, **k):
    if len(a) >= 2 and np.ndim(a[0]) == 1 and len(a[0]) > 5:
        REC.append(dict(x=list(map(float, a[0])), y=list(map(float, a[1])), label=k.get('label')))
    return _p(*a, **k)
plt.plot = p
atexit.register(lambda: json.dump(REC, open('scan_capture.json', 'w'), indent=1))
sys.path.insert(0, D); sys.argv = [D + '/doFit.py'] + sys.argv[1:]
runpy.run_path(D + '/doFit.py', run_name='__main__')
