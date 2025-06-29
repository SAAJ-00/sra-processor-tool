import os
from pathlib import Path

DEFAULT_CONFIG = {
    'threads': 8,
    'max_size': '30G',
    'output_dir': Path('.').absolute(),
    'keep_temp': False,
    'trim_params': {
        # Parámetros para lecturas cortas
        'quality_phred': 30,
        'min_length': 50,
        'cut_window_size': 4,
        'cut_mean_quality': 25,
        'disable_adapter_trimming': False,
        'disable_quality_filtering': False,
        'disable_length_filtering': False,

        # Parámetros específicos para lecturas largas
        'long_read_settings': {
            'min_quality': 10,
            'min_length': 1000,
            'disable_adapter_trimming': False,
            'disable_quality_filtering': False
        }
    }
}