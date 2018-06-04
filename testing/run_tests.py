from __future__ import print_function
import glob
import sys
import genomeweb as gw

files = [f for f in glob.glob('genomes/*.fna')]
if len(sys.argv) > 1:
    if sys.argv[1] == '-v':
        quiet = False
    else:
        print('Unkown argument: %s' % sys.argv[1])
        print('Use "-v" for verbose output.')
        sys.exit(1)
else:
    quiet = True

def test(func, i):
    pf = 'Fail'
    e = ''
    print('Running Test %d...' % i)
    try:
        func()
        pf = 'Pass'
    except Exception as ex:
        e = ex
    print('Test %d:\t%s\t%s' % (i, pf, e))
    return e

opts = dict(
    working_directory='scratch',
    matches_opts=dict(
        quiet=quiet),
    reorder_opts=dict(
        quiet=quiet),
    quiet=quiet
    )
# TODO: update tests
tests = [
    lambda: gw.create_web(files[1:], reference_genome=files[0], **opts),
    lambda: gw.create_web(files[1:], reference_genome=[
        files[0], '', ''], **opts)]
fails = []
for i, t in enumerate(tests):
    e = test(t, i + 1)
    if e:
        fails.append(i + 1)
if fails:
    print('Not all tests passed, check install.')
    sys.exit(1)
else:
    print('All tests passed.')
    sys.exit(0)
