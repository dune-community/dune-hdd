#!/usr/bin/env python

from generic import dune_module, examples

Example = examples[3]['aluconformgrid']['fem']['istl']
spe10_cfg = {'filename': 'spe_perm.dat',
             'upper_right': '[2 5 1]'}

logger_cfg = Example.logger_options()
logger_cfg.set('info_color', 'blue', True)
logger_cfg.set('info', 99, True)

grid_cfg = Example.grid_options('stuff.grid.provider.cube')
grid_cfg.set('type', 'stuff.grid.provider.cube', True)
grid_cfg.set('num_elements', '[12 44 16]', True)
grid_cfg.set('upper_right', spe10_cfg['upper_right'], True)

boundary_cfg = dune_module.Dune.Stuff.Common.Configuration()
boundary_cfg.set('type', 'stuff.grid.boundaryinfo.normalbased')
boundary_cfg.set('default', 'neumann')
boundary_cfg.set('dirichlet.0', '[0 1 0]')

problem_cfg = dune_module.Dune.Stuff.Common.Configuration()
problem_cfg.set('type', 'hdd.linearelliptic.problem.default')
problem_cfg.set('diffusion_factor.type', 'stuff.function.constant')
problem_cfg.set('diffusion_factor.name', 'diffusion_factor')
problem_cfg.set('diffusion_factor.value', '1')
problem_cfg.set('diffusion_tensor.type', 'pymor.function.spe10.model2')
problem_cfg.set('diffusion_tensor.upper_right', spe10_cfg['upper_right'])
problem_cfg.set('diffusion_tensor.blockade_width', '0.2')
problem_cfg.set('diffusion_tensor.blockade_value', '1')
problem_cfg.set('force.type', 'pymor.function.affinelydecomposabledefault')
problem_cfg.set('force.name', 'force')
problem_cfg.set('force.component.0.type', 'stuff.function.indicator')
problem_cfg.set('force.component.0.name', 'sink')
problem_cfg.set('force.component.0.0.domain', '[0.167 0.5; 1.82 2.05; 0 1]')
problem_cfg.set('force.component.0.0.value', '-1000')
problem_cfg.set('force.coefficient.0.sink', '1')
problem_cfg.set('force.coefficient.0.expression', 'sink[0]')
problem_cfg.set('dirichlet.type', 'stuff.function.constant')
problem_cfg.set('dirichlet.name', 'dirichlet')
problem_cfg.set('dirichlet.value', '0')
problem_cfg.set('neumann.type', 'stuff.function.indicator')
problem_cfg.set('neumann.name', 'neumann')
problem_cfg.set('neumann.0.domain', '[-999 999; -999 0.03; -999 999]')
problem_cfg.set('neumann.0.value', '1')

example = Example(logger_cfg, grid_cfg, boundary_cfg, problem_cfg)
# example.visualize_grid('spe10.grid')

solver_tp = 'bicgstab.amg.ilu0'
solver_cfg = solver_cfg = Example.solver_options(solver_tp)
solver_cfg.set('preconditioner.anisotropy_dim', 3, True)
solver_cfg.set('preconditioner.isotropy_dim', 3, True)

d = example.discretization()

for mu_ in ((0, 0), (0, 1), (1, 0), (1, 1)):
    mu = dune_module.Dune.Pymor.Parameter(d.parameter_type(),
                                          list((list((mu_[0],)), list((mu_[1],)))))
    # example.visualize_problem('spe10.problem.' + mu.report_for_filename(), mu)
    U = d.create_vector()
    d.solve(solver_cfg, U, mu)
    d.visualize(U,
                'spe10.pressure.' + mu.report_for_filename(),
                'pressure' + mu.report_for_filename())
    # example.visualize_darcy_velocity(U, 'spe10.velocity.' + mu.report_for_filename(), 'velocity', mu)

