import os, sys, shutil
import subprocess as s

filename_lm     = 'models/ball_and_stick.lm'
filename_run    = 'results/run_ball_and_stick.lm'

if os.path.isfile(filename_run):
    os.remove(filename_run)
    print('Run file was removed.')

print('Create run file: ', filename_run)
os.makedirs('results', exist_ok=True)
shutil.copy(filename_lm, filename_run)

com = ['lm','-r', '1', '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename_run]
print(' '.join(com))
s.call(com)