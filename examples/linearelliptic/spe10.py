
from generic import dune_module, examples

Example = examples[3]['yaspgrid']['pdelab']['istl']

spe10_cfg = {'filename': 'spe_perm.dat',
             'upper_right': '[365.76 670.56 51.816]',
             'num_elements': '[60 220 85]'}

logger_cfg = Example.logger_options()
logger_cfg.set('info_color', 'blue', True)
logger_cfg.set('info', 99, True)

grid_cfg = Example.grid_options('stuff.grid.provider.cube')
grid_cfg.set('type', 'stuff.grid.provider.cube', True)
grid_cfg.set('num_elements', spe10_cfg['num_elements'], True)
grid_cfg.set('upper_right', spe10_cfg['upper_right'], True)

boundary_cfg = Example.boundary_options('stuff.grid.boundaryinfo.alldirichlet')

problem_cfg = dune_module.Dune.Stuff.Common.Configuration()
problem_cfg.set('type', 'hdd.linearelliptic.problem.default')
problem_cfg.set('diffusion_factor.type', 'stuff.function.constant')
problem_cfg.set('diffusion_factor.name', 'diffusion_factor')
problem_cfg.set('diffusion_factor.value', '1')
problem_cfg.set('diffusion_tensor.type', 'stuff.function.spe10.model2')
problem_cfg.set('force.type', 'stuff.function.constant')
problem_cfg.set('force.name', 'force')
problem_cfg.set('force.value', '1')
problem_cfg.set('dirichlet.type', 'stuff.function.constant')
problem_cfg.set('dirichlet.name', 'dirichlet')
problem_cfg.set('dirichlet.value', '0')
problem_cfg.set('neumann.type', 'stuff.function.constant')
problem_cfg.set('neumann.name', 'neumann')
problem_cfg.set('neumann.value', '0')

example = Example(logger_cfg, grid_cfg, boundary_cfg, problem_cfg)
example.visualize('spe10_example')

