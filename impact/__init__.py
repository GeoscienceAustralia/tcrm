"""
Hazimp integration for TCRM


Send the windfield hazard to hazimp, which will combine it with
exposure and vulnerabilty data to report on resulting impact of the cyclone.
"""

# Minimalist integration approach:
#   Assume hazimp installed, invoke from command line
#   (rather than importing API).

import subprocess
import tempfile
import logging
import yaml
import glob
import os.path

def impact(gust, exposure, output='impact.csv'):
    """Invoke hazimp for a particular gust file"""

    settings = [{'template': 'wind_v3'},
                {'load_exposure': dict(file_name=exposure,
                                       exposure_latitude='LATITUDE',
                                       exposure_longitude='LONGITUDE')},
                {'load_wind_ascii': 'NetCDF:{}:vmax'.format(gust)},
                {'calc_struct_loss':
                    dict(replacement_value_label='REPLACEMENT_VALUE')},
                {'save': output}]
    hazimp_config = yaml.dump(settings)
    logging.debug(hazimp_config)

    with tempfile.NamedTemporaryFile(mode='w') as f:
        f.write(hazimp_config)
        f.flush()
        #f.close()

        cmd = ['hazimp', '-c', str(f.name)]
        logging.debug(' '.join(cmd))

        status = subprocess.call(cmd)
        if status:
            raise Exception("hazimp failure: {}".format(settings))

def run(workDir, exposure):
    """Calculate impact for each available gust"""

    inputpattern = os.path.join(workDir, 'windfield', 'gust.*.nc')
    logging.info('Preparing to calculate impact for {}'.format(inputpattern))

    # Confirm output location exists:
    #os.makedirs(os.path.join(workDir, 'impact'), exist_ok=True) #py3
    try:
        os.makedirs(os.path.join(workDir, 'impact'))
    except OSError as e:
        import errno
        if e.errno != errno.EEXIST:
            raise

    # Could balance the loop here over parallel resources,
    # or could leave hazimp to implement any parallelisation.
    for gust in glob.glob(inputpattern):
        output = os.path.join(workDir, 'impact',
                              os.path.basename(gust)[:-3] + '.csv')
        impact(gust, exposure, output)

def run_optional(config):
    """Interpret config and (if so instructed) then calculate impact"""

    if config.has_section('Actions'):
        action = config.getboolean('Actions', 'ExecuteImpact', fallback=False)
    else:
        action = config.has_section('Impact')

    if action:
        run(config.get('Output', 'Path'), config.get('Impact', 'Exposure'))
