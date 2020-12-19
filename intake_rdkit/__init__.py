#-----------------------------------------------------------------------------
# Copyright (c) 2020, Greg Landrum
# All rights reserved.
#
# The full license is in the LICENSE file, distributed with this software.
#-----------------------------------------------------------------------------
__version__ = "0.1.0"

import intake  # Import this first to avoid circular imports during discovery.
from intake.container import register_container
from .sdf import SDFSource

