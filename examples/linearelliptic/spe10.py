#!/usr/bin/env python

from generic import dune_module, examples

DuneException = dune_module.Dune.Exception

for la_backend, Example in examples[3]['aluconformgrid']['fem'].items():
    spe10_cfg = {'filename': 'spe_perm.dat',
                 'upper_right': '[2 5 1]'}

    logger_cfg = Example.logger_options()
    logger_cfg.set('info_color', 'blue', True)
    logger_cfg.set('info', 99, True)

    grid_cfg = Example.grid_options('stuff.grid.provider.cube')
    grid_cfg.set('type', 'stuff.grid.provider.cube', True)
    grid_cfg.set('num_elements', '[16 40 8]', True)
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
    problem_cfg.set('diffusion_tensor.type', 'stuff.function.spe10.model2')
    problem_cfg.set('diffusion_tensor.upper_right', spe10_cfg['upper_right'])
    problem_cfg.set('force.type', 'stuff.function.constant')
    problem_cfg.set('force.name', 'force')
    problem_cfg.set('force.value', '0')
    problem_cfg.set('dirichlet.type', 'stuff.function.constant')
    problem_cfg.set('dirichlet.name', 'dirichlet')
    problem_cfg.set('dirichlet.value', '0')
    problem_cfg.set('neumann.type', 'stuff.function.indicator')
    problem_cfg.set('neumann.name', 'neumann')
    problem_cfg.set('neumann.0.domain', '[-999 999; -999 0.1; -999 999]')
    problem_cfg.set('neumann.0.value', '1')

    example = Example(logger_cfg, grid_cfg, boundary_cfg, problem_cfg)
    example.visualize('spe10_example')

    for solver_tp in list(Example.solver_options()):
        if solver_tp in ('bicgstab.ilut', 'lu.sparse', 'superlu', 'bicgstab.amg.ilu0'):
            solver_cfg = solver_cfg = Example.solver_options(solver_tp)
            solver_cfg.set('preconditioner.anisotropy_dim', 3, True)
            solver_cfg.set('preconditioner.isotropy_dim', 3, True)

            cg_d = example.cg_discretization()
            cg_U = cg_d.create_vector()
            try:
                cg_d.solve(solver_cfg, cg_U, dune_module.Dune.Pymor.Parameter())
                solver_cfg.report()
                prefix = ''
            except DuneException as ee:
                prefix = 'failed__'
                print(ee)
            cg_d.visualize(cg_U, 'spe10_example.solution.cg.{}{}{}'.format(la_backend, prefix, solver_tp), 'pressure_cg')

            swipdg_d = example.swipdg_discretization()
            swipdg_U = swipdg_d.create_vector()
            try:
                swipdg_d.solve(solver_cfg, swipdg_U, dune_module.Dune.Pymor.Parameter())
                solver_cfg.report()
                prefix = ''
            except DuneException as ee:
                prefix = 'failed__'
                print(ee)
            swipdg_d.visualize(swipdg_U, 'spe10_example.solution.swipdg.{}{}{}'.format(la_backend, prefix, solver_tp), 'pressure_dg')
