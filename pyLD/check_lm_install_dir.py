import os, sys, site
from pyLD import locate_cuda
VALID_CUDA_VERSIONS = ['11.0','11.4']

def binaries_directory(user):
	"""Return the installation directory, or None"""
	if user == True:
		paths = (site.getusersitepackages(),)
	else:
		py_version = '%s.%s' % (sys.version_info[0], sys.version_info[1])
		paths = (s % (py_version) for s in (
			sys.prefix + '/lib/python%s/dist-packages/',
			sys.prefix + '/lib/python%s/site-packages/',
			sys.prefix + '/local/lib/python%s/dist-packages/',
			sys.prefix + '/local/lib/python%s/site-packages/',
			'/Library/Python/%s/site-packages/',
		))

	valid_path_list = []
	for path in paths:
		if os.path.exists(path):
			valid_path_list.append(path)
	return valid_path_list


def check_lm_install_dir():
	if os.name== 'posix':
		cuda_version = locate_cuda()
		print('Detected CUDA version: ', cuda_version)
	else:
		print('Linux CUDA undetected.')
	print('Valid CUDA versions : ', ','.join(VALID_CUDA_VERSIONS) )

	module_paths1 = binaries_directory(True)
	module_paths2 = binaries_directory(False)
	module_paths  = module_paths1 + module_paths2
	for module_path in module_paths:
		bin_dir = os.path.join(module_path,'lm','bin')
		lib_dir = os.path.join(module_path,'lm','lib')
		if os.path.exists(bin_dir) and os.path.exists(lib_dir):
			print('Installed $LM_BIN_PATH: ', bin_dir)
			print('Installed $LM_LIB_PATH: ', lib_dir)
			sys.exit()
	print('The pair of lm/bin and lm/lib was not found.')


